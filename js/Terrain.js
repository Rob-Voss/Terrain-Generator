"use strict";

var rNorm = (function () {
  var z2 = null;

  function rnorm() {
    if (z2 != null) {
      var tmp = z2;
      z2 = null;
      return tmp;
    }
    var x1 = 0,
      x2 = 0,
      w = 2.0;
    while (w >= 1) {
      x1 = runIf(-1, 1);
      x2 = runIf(-1, 1);
      w = x1 * x1 + x2 * x2;
    }
    w = Math.sqrt(-2 * Math.log(w) / w);
    z2 = x2 * w;

    return x1 * w;
  }

  return rnorm;
})();

var city_counts = {
  'shore': (12, 20),
  'island': (5, 10),
  'mountain': (10, 25),
  'desert': (5, 10)
},
terr_counts = {
  'shore': (3, 7),
  'island': (2, 4),
  'mountain': (3, 6),
  'desert': (3, 5)
},
riverpercs = {
  'shore': 5,
  'island': 3,
  'mountain': 8,
  'desert': 1
};

const
  /**
   * Default font sizes for the map
   *
   * @typedef {Object} fontSizeObject
   * @property {number} region
   * @property {number} city
   * @property {number} town
   */
  fontSizeObject = {
    region: 40,
    city: 25,
    town: 20
  },

  /**
   * Extent Object
   *
   * @typedef {Object} extentObject
   * @property {number} width
   * @property {number} height
   */
  extentObject = {
    width: 1,
    height: 1
  },

  /**
   * Default options for the map generator
   *
   * @typedef {Object} paramsObject
   * @property {extentObject} extent
   * @property {function} generator
   * @property {number} width
   * @property {number} height
   * @property {number} numPts
   * @property {number} numCities
   * @property {number} numTerritories
   * @property {fontSizeObject} fontSizes
   */
  paramsObject = {
    extent: extentObject,
    generator: generateCoastLine,
    fontSizes: fontSizeObject,
    width: 1000,
    height: 1000,
    numPts: 10000,
    numCities: 15,
    numTerritories: 5
  },

  /**
   * Mesh object
   *
   * @typeDef {Object} meshObject
   * @property {Array} adj
   * @property {Array} edges
   * @property {extentObject} extent
   * @property {function} map
   * @property {Array} pts
   * @property {Array} tris
   * @property {d3.voronoi} vor
   * @property {Array} vxs
   */
  meshObject = {
    adj: [],
    edges: [],
    extent: {},
    map: {},
    pts: [],
    tris: [],
    vor: {},
    vxs: []
  };

/**
 *
 * @param {paramsObject} params
 */
function generateCoastLine(params = paramsObject) {
  var mesh = Mesh.generateGoodMesh(params.numPts, params.extent);
  var landscape = add(
    Mesh.slope(mesh, randomVector(4)),
    Mesh.cone(mesh, runIf(-1, -1)),
    Mesh.mountains(mesh, 150, 0.2)
  );

  for (let i = 0; i < 10; i++) {
    landscape = Mesh.relax(landscape);
  }

  landscape = Mesh.doErosion(landscape, runIf(0, 0.1), 5);
  landscape = Mesh.setSeaLevel(landscape, runIf(0.2, 0.6));
  landscape = Mesh.fillSinks(landscape);
  landscape = Mesh.cleanCoast(landscape, 3);

  return landscape;
}

/**
 *
 * @param params
 */
function generateLake(params = paramsObject) {
  var mesh = Mesh.generateGoodMesh(params.numPts, params.extent);

  var landscape = add(
    Mesh.cone(mesh, 2),
    Mesh.mountains(mesh, 50)
  );

  for (let i = 0; i < 10; i++) {
    landscape = Mesh.relax(landscape);
  }

  landscape = Mesh.doErosion(landscape, runIf(0, 0.1), 5);
  landscape = Mesh.setSeaLevel(landscape, params.seaLevel || 0.5);
  landscape = Mesh.fillSinks(landscape);
  landscape = Mesh.cleanCoast(landscape, 3);

  return landscape;
}

/**
 *
 * @param params
 */
function generateIsland(params = paramsObject) {
  var mesh = Mesh.generateGoodMesh(params.numPts, params.extent);

  var landscape = add(
    Mesh.cone(mesh, -1),
    Mesh.mountains(mesh, 50)
  );

  for (let i = 0; i < 10; i++) {
    landscape = Mesh.relax(landscape);
  }

  landscape = Mesh.peaky(landscape);
  landscape = Mesh.doErosion(landscape, runIf(0, 0.1), 5);
  landscape = Mesh.setSeaLevel(landscape, params.seaLevel || 0.5);
  landscape = Mesh.fillSinks(landscape);
  landscape = Mesh.cleanCoast(landscape, 3);

  return landscape;
}



function randomVector(scale) {
  return [scale * rNorm(), scale * rNorm()];
}

function runIf(lo, hi) {
  return lo + Math.random() * (hi - lo);
}

function add() {
  var n = arguments[0].length;
  var newVals = Mesh.zero(arguments[0].mesh);
  for (var i = 0; i < n; i++) {
    for (var j = 0; j < arguments.length; j++) {
      newVals[i] += arguments[j][i];
    }
  }

  return newVals;
}



/**
 *
 */
class Mesh {

  /**
   * The point where the three medians of the triangle intersect
   *
   * @param {array} pts
   * @returns {*[]}
   */
  static centroid(pts) {
    var x = 0, y = 0;
    for (let i = 0; i < pts.length; i++) {
      x += pts[i][0];
      y += pts[i][1];
    }

    return [x / pts.length, y / pts.length];
  }

  /**
   * Smooth out the lines of the coast path points by repeatedly applying
   * a filter where points which are below sea level, but a majority of
   * whose neighbors are above sea level, get pulled up, and vice versa
   * for points which are above sea level and have undersea neighbors.
   *
   * @param {object} landscape
   * @param {number} iterations
   * @returns {*}
   */
  static cleanCoast(landscape, iterations) {
    for (let it = 0; it < iterations; it++) {
      var changed = 0,
        newH = Mesh.zero(landscape.mesh);
      for (let i = 0; i < landscape.length; i++) {
        newH[i] = landscape[i];
        let neighbors = Mesh.neighbors(landscape.mesh, i),
          count = 0,
          best = -999999;
        if (landscape[i] <= 0 || neighbors.length != 3) {
          continue;
        }

        for (let j = 0; j < neighbors.length; j++) {
          if (landscape[neighbors[j]] > 0) {
            count++;
          } else if (landscape[neighbors[j]] > best) {
            best = landscape[neighbors[j]];
          }
        }
        if (count > 1) {
          continue;
        }
        newH[i] = best / 2;
        changed++;
      }
      landscape = newH;
      newH = Mesh.zero(landscape.mesh);
      for (let i = 0; i < landscape.length; i++) {
        newH[i] = landscape[i];
        let nbs = Mesh.neighbors(landscape.mesh, i),
          count = 0,
          best = 999999;
        if (landscape[i] > 0 || nbs.length != 3) {
          continue;
        }

        for (let j = 0; j < nbs.length; j++) {
          if (landscape[nbs[j]] <= 0) {
            count++;
          } else if (landscape[nbs[j]] < best) {
            best = landscape[nbs[j]];
          }
        }
        if (count > 1) {
          continue;
        }
        newH[i] = best / 2;
        changed++;
      }
      landscape = newH;
    }

    return landscape;
  }

  /**
   * Generate a cone or an inverted cone to make a hole
   *
   * @param {meshObject} mesh
   * @param {number} slope
   */
  static cone(mesh, slope) {
    return mesh.map(function (x) {
      return Math.pow(x[0] * x[0] + x[1] * x[1], 0.5) * slope;
    });
  }

  /**
   * Generate a contoured path at a specific level
   *
   * @param {object} landscape
   * @param {number} level
   * @returns {Array}
   */
  static contour(landscape, level = 0) {
    var edges = [];
    for (let i = 0; i < landscape.mesh.edges.length; i++) {
      var e = landscape.mesh.edges[i],
        e0 = e[0],
        e1 = e[1],
        lE0 = landscape[e0],
        lE1 = landscape[e1];
      if (e[3] == undefined) {
        continue;
      }
      if (Mesh.isNearEdge(landscape.mesh, e0) || Mesh.isNearEdge(landscape.mesh, e1)) {
        continue;
      }
      if ((lE0 > level && lE1 <= level) ||
        (lE1 > level && lE0 <= level)) {
        edges.push([e[2], e[3]]);
      }
    }

    return Mesh.mergeSegments(edges);
  }

  /**
   * Get the distance between two points in the mesh?
   *
   * @param {meshObject} mesh
   * @param {number} i
   * @param {number} j
   * @returns {number}
   */
  static distance(mesh, i, j) {
    var p = mesh.vxs[i],
      q = mesh.vxs[j];

    return Math.sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]));
  }

  /**
   * Erode the coast line of the landscape
   *
   * @param {object} landscape
   * @param {number} amount
   * @param {number} n
   * @returns {object} landscape
   */
  static doErosion(landscape, amount = 0.1, n = 1) {
    landscape = Mesh.fillSinks(landscape);
    for (let i = 0; i < n; i++) {
      landscape = Mesh.erode(landscape, amount);
      landscape = Mesh.fillSinks(landscape);
    }

    return landscape;
  }

  /**
   * Generate the hills for the rivers?
   *
   * @param {object} landscape
   * @returns {Array} downs
   */
  static downhill(landscape) {
    if (landscape.downhill) {
      return landscape.downhill;
    }

    function downFrom(i) {
      if (Mesh.isEdge(landscape.mesh, i)) {
        return -2;
      }
      var best = -1,
        bestH = landscape[i],
        nbs = Mesh.neighbors(landscape.mesh, i);
      for (let j = 0; j < nbs.length; j++) {
        if (landscape[nbs[j]] < bestH) {
          bestH = landscape[nbs[j]];
          best = nbs[j];
        }
      }

      return best;
    }

    var downs = [];
    for (let i = 0; i < landscape.length; i++) {
      downs[i] = downFrom(i);
    }
    landscape.downhill = downs;

    return downs;
  }

  /**
   *
   * @param {object} landscape
   * @param {number} p
   */
  static dropEdge(landscape, p) {
    p = p || 4;
    var newH = Mesh.zero(landscape.mesh);
    for (let i = 0; i < landscape.length; i++) {
      var v = landscape.mesh.vxs[i],
        x = 2.4 * v[0] / landscape.mesh.extent.width,
        y = 2.4 * v[1] / landscape.mesh.extent.height;
      newH[i] = landscape[i] - Math.exp(10 * (Math.pow(Math.pow(x, p) + Math.pow(y, p), 1 / p) - 1));
    }

    return newH;
  }

  /**
   * Generate a new landscape with eroded coast lines?
   *
   * @param {object} landscape
   * @param {number} amount
   * @returns {meshObject} newH
   */
  static erode(landscape, amount) {
    var er = Mesh.erosionRate(landscape),
      newH = Mesh.zero(landscape.mesh),
      maxR = d3.max(er);
    for (let i = 0; i < landscape.length; i++) {
      newH[i] = landscape[i] - amount * (er[i] / maxR);
    }

    return newH;
  }

  /**
   * Set and erode the landscape rate?
   *
   * @param {object} landscape
   * @returns {meshObject} newH
   */
  static erosionRate(landscape) {
    var flux = Mesh.getFlux(landscape),
      slope = Mesh.getSlope(landscape),
      newH = Mesh.zero(landscape.mesh);
    for (let i = 0; i < landscape.length; i++) {
      var river = Math.sqrt(flux[i]) * slope[i],
        creep = slope[i] * slope[i],
        total = 1000 * river + creep;
      total = total > 200 ? 200 : total;
      newH[i] = total;
    }

    return newH;
  }

  /**
   * Planchon-Darboux algorithm
   *
   * @param {object} landscape
   * @param {number} epsilon
   * @returns {meshObject} newH
   */
  static fillSinks(landscape, epsilon = 1e-5) {
    var infinity = 999999,
      newH = Mesh.zero(landscape.mesh);
    for (let i = 0; i < landscape.length; i++) {
      if (Mesh.isNearEdge(landscape.mesh, i)) {
        newH[i] = landscape[i];
      } else {
        newH[i] = infinity;
      }
    }

    while (true) {
      var changed = false;
      for (let i = 0; i < landscape.length; i++) {
        if (newH[i] == landscape[i]) {
          continue;
        }
        var nbs = Mesh.neighbors(landscape.mesh, i);
        for (let j = 0; j < nbs.length; j++) {
          if (landscape[i] >= newH[nbs[j]] + epsilon) {
            newH[i] = landscape[i];
            changed = true;
            break;
          }
          var oh = newH[nbs[j]] + epsilon;
          if ((newH[i] > oh) && (oh > landscape[i])) {
            newH[i] = oh;
            changed = true;
          }
        }
      }
      if (!changed) {
        return newH;
      }
    }
  }

  /**
   *
   * @param {object} landscape
   */
  static findSinks(landscape) {
    var dh = Mesh.downhill(landscape),
      sinks = [];
    for (let i = 0; i < dh.length; i++) {
      var node = i;
      while (true) {
        if (Mesh.isEdge(landscape.mesh, node)) {
          sinks[i] = -2;
          break;
        }
        if (dh[node] == -1) {
          sinks[i] = node;
          break;
        }
        node = dh[node];
      }
    }
  }

  /**
   *
   * @param {number} n
   * @param {extentObject} extent
   * @returns {meshObject}
   */
  static generateGoodMesh(n, extent) {
    n = n || 1;
    extent = extent || extentObject;
    var pts = Mesh.generateGoodPoints(n, extent);

    return Mesh.makeMesh(pts, extent);
  }

  /**
   *
   * @param {number} n
   * @param {extentObject} extent
   * @returns {Array}
   */
  static generateGoodPoints(n, extent) {
    n = n || 1;
    extent = extent || extentObject;
    var pts = Mesh.generatePoints(n, extent);
    pts = pts.sort(function (a, b) {
      return a[0] - b[0];
    });

    return Mesh.improvePoints(pts, 1, extent);
  }

  /**
   *
   * @param {number} n
   * @param {extentObject} extent
   * @returns {Array} pts
   */
  static generatePoints(n, extent) {
    n = n || 1;
    extent = extent || extentObject;
    var pts = [];
    for (let i = 0; i < n; i++) {
      pts.push([
        (Math.random() - 0.5) * extent.width,
        (Math.random() - 0.5) * extent.height
      ]);
    }

    return pts;
  }

  /**
   * Get the water flux from each point to its 'downhill point'.
   *
   * @param {object} landscape
   * @returns {meshObject} flux
   */
  static getFlux(landscape) {
    var dh = Mesh.downhill(landscape),
      idxs = [],
      flux = Mesh.zero(landscape.mesh);

    for (let i = 0; i < landscape.length; i++) {
      idxs[i] = i;
      flux[i] = 1 / landscape.length;
    }

    idxs.sort(function (a, b) {
      return landscape[b] - landscape[a];
    });

    for (let i = 0; i < landscape.length; i++) {
      var j = idxs[i];
      if (dh[j] >= 0) {
        flux[dh[j]] += flux[j];
      }
    }

    return flux;
  }

  /**
   * Get the slope values
   *
   * @param {object} landscape
   * @returns {meshObject} slope
   */
  static getSlope(landscape) {
    var dh = Mesh.downhill(landscape),
      slope = Mesh.zero(landscape.mesh);
    for (let i = 0; i < landscape.length; i++) {
      var s = Mesh.triSlope(landscape, i);
      slope[i] = Math.sqrt(s[0] * s[0] + s[1] * s[1]);
      // continue;
      if (dh[i] < 0) {
        slope[i] = 0;
      } else {
        slope[i] = (landscape[i] - landscape[dh[i]]) / Mesh.distance(landscape.mesh, i, dh[i]);
      }
    }
    return slope;
  }

  /**
   * Lloyd relaxation to improve the point set
   * @param {Array} pts
   * @param {number} n
   * @param {extentObject} extent
   * @returns {Array} pts
   */
  static improvePoints(pts, n, extent) {
    n = n || 1;
    extent = extent || extentObject;
    for (let i = 0; i < n; i++) {
      pts = Mesh.voronoi(pts, extent).polygons(pts).map(Mesh.centroid);
    }

    return pts;
  }

  /**
   * Check if it is an edge
   * @param {meshObject} mesh
   * @param {number} i
   * @returns {boolean}
   */
  static isEdge(mesh, i) {
    return (mesh.adj[i].length < 3);
  }

  /**
   * Check if it is near an edge
   * @param {meshObject} mesh
   * @param {number} i
   * @returns {boolean}
   */
  static isNearEdge(mesh, i) {
    var x = mesh.vxs[i][0],
      y = mesh.vxs[i][1],
      w = mesh.extent.width,
      h = mesh.extent.height;

    return x < -0.45 * w || x > 0.45 * w || y < -0.45 * h || y > 0.45 * h;
  }

  /**
   *
   * @param {Array} points
   * @param {extentObject} extent
   * @returns {meshObject} mesh
   */
  static makeMesh(points, extent) {
    var voronoi = Mesh.voronoi(points, extent),
      vxs = [],
      vxIds = {},
      adjacencies = [],
      edges = [],
      triangles = [];
    for (let i = 0; i < voronoi.edges.length; i++) {
      var edge = voronoi.edges[i];
      if (edge == undefined) {
        continue;
      }
      var e0 = vxIds[edge[0]],
        e1 = vxIds[edge[1]];
      if (e0 == undefined) {
        e0 = vxs.length;
        vxIds[edge[0]] = e0;
        vxs.push(edge[0]);
      }
      if (e1 == undefined) {
        e1 = vxs.length;
        vxIds[edge[1]] = e1;
        vxs.push(edge[1]);
      }
      adjacencies[e0] = adjacencies[e0] || [];
      adjacencies[e0].push(e1);
      adjacencies[e1] = adjacencies[e1] || [];
      adjacencies[e1].push(e0);
      edges.push([e0, e1, edge.left, edge.right]);
      triangles[e0] = triangles[e0] || [];
      if (!triangles[e0].includes(edge.left)) {
        triangles[e0].push(edge.left);
      }
      if (edge.right && !triangles[e0].includes(edge.right)) {
        triangles[e0].push(edge.right);
      }
      triangles[e1] = triangles[e1] || [];
      if (!triangles[e1].includes(edge.left)) {
        triangles[e1].push(edge.left);
      }
      if (edge.right && !triangles[e1].includes(edge.right)) {
        triangles[e1].push(edge.right);
      }
    }

    /**
     * @type {meshObject} mesh
     */
    var mesh = {
      adj: adjacencies,
      edges: edges,
      extent: extent,
      pts: points,
      tris: triangles,
      vor: voronoi,
      vxs: vxs
    };

    mesh.map = function (f) {
      var mapped = vxs.map(f);
      mapped.mesh = mesh;

      return mapped;
    };

    return mesh;
  }

  /**
   *
   * @param {object} landscape
   * @param {function} f
   * @returns {object} newH
   */
  static map(landscape, f) {
    var newH = landscape.map(f);
    newH.mesh = landscape.mesh;

    return newH;
  }

  /**
   *
   * @param {Array} segments
   * @returns {Array} paths
   */
  static mergeSegments(segments) {
    var adj = {};
    for (let i = 0; i < segments.length; i++) {
      var seg = segments[i],
        a0 = adj[seg[0]] || [],
        a1 = adj[seg[1]] || [];
      a0.push(seg[1]);
      a1.push(seg[0]);
      adj[seg[0]] = a0;
      adj[seg[1]] = a1;
    }
    var done = [],
      paths = [],
      path = null;
    while (true) {
      if (path === null) {
        for (let i = 0; i < segments.length; i++) {
          if (done[i]) continue;
          done[i] = true;
          path = [segments[i][0], segments[i][1]];
          break;
        }
        if (path === null) {
          break;
        }
      }
      var changed = false;
      for (let i = 0; i < segments.length; i++) {
        if (done[i]) {
          continue;
        }
        if (adj[path[0]].length === 2 && segments[i][0] === path[0]) {
          path.unshift(segments[i][1]);
        } else if (adj[path[0]].length === 2 && segments[i][1] === path[0]) {
          path.unshift(segments[i][0]);
        } else if (adj[path[path.length - 1]].length === 2 && segments[i][0] === path[path.length - 1]) {
          path.push(segments[i][1]);
        } else if (adj[path[path.length - 1]].length === 2 && segments[i][1] === path[path.length - 1]) {
          path.push(segments[i][0]);
        } else {
          continue;
        }
        done[i] = true;
        changed = true;
        break;
      }
      if (!changed) {
        paths.push(path);
        path = null;
      }
    }

    return paths;
  }

  /**
   * Generate some mountains
   *
   * @param {meshObject} mesh
   * @param {number} number
   * @param {number} radius
   * @returns {meshObject} newVals
   */
  static mountains(mesh, number, radius = 0.05) {
    let mounts = [],
      newVals = Mesh.zero(mesh);
    for (let i = 0; i < number; i++) {
      let x = mesh.extent.width * (Math.random() - 0.5),
        y = mesh.extent.height * (Math.random() - 0.5);
      mounts.push([x, y]);
    }

    for (let i = 0; i < mesh.vxs.length; i++) {
      let p = mesh.vxs[i];
      for (let j = 0; j < number; j++) {
        let m = mounts[j],
          val = Math.pow(Math.exp(-((p[0] - m[0]) * (p[0] - m[0]) + (p[1] - m[1]) * (p[1] - m[1])) / (2 * radius * radius)), 2);
        newVals[i] += val;
      }
    }

    return newVals;
  }

  /**
   *
   * @param {meshObject} mesh
   * @param {number} i
   * @returns {Array} nbs
   */
  static neighbors(mesh, i) {
    var onbs = mesh.adj[i],
      nbs = [];
    for (let i = 0; i < onbs.length; i++) {
      nbs.push(onbs[i]);
    }

    return nbs;
  }

  /**
   * Rescale the heights to lie in the range 0-1
   *
   * @param {object} landscape
   * @returns {*}
   */
  static normalize(landscape) {
    var lo = d3.min(landscape),
      hi = d3.max(landscape);

    return Mesh.map(landscape, function (x) {
      return (x - lo) / (hi - lo);
    });
  }

  /**
   * Normalize, then take the square root of the height value
   * to round off the tops of hills
   *
   * @param {object} landscape
   * @returns {*}
   */
  static peaky(landscape) {
    return Mesh.map(Mesh.normalize(landscape), Math.sqrt);
  }

  /**
   *
   * @param {object} landscape
   * @param {number} q
   * @returns {number}
   */
  static quantile(landscape, q) {
    var sortedH = [];
    for (let i = 0; i < landscape.length; i++) {
      sortedH[i] = landscape[i];
    }
    sortedH.sort(d3.ascending);

    return d3.quantile(sortedH, q);
  }

  /**
   * Replace each height value with the average of its neighbours
   * to smooth the surface
   *
   * @param {object} landscape
   * @returns {meshObject} newH
   */
  static relax(landscape) {
    var newH = Mesh.zero(landscape.mesh);
    for (let i = 0; i < landscape.length; i++) {
      var nbs = Mesh.neighbors(landscape.mesh, i);
      if (nbs.length < 3) {
        newH[i] = 0;
        continue;
      }
      newH[i] = d3.mean(nbs.map(function (j) {
        return landscape[j];
      }));
    }

    return newH;
  }

  /**
   * Relax the points in the middle of the path towards their upstream
   * and downstream neighbors
   *
   * @param {Array} path
   * @returns {Array} newPath
   */
  static relaxPath(path) {
    var newPath = [path[0]];
    for (let i = 1; i < path.length - 1; i++) {
      var newPt = [
        0.25 * path[i - 1][0] + 0.5 * path[i][0] + 0.25 * path[i + 1][0],
        0.25 * path[i - 1][1] + 0.5 * path[i][1] + 0.25 * path[i + 1][1]
      ];
      newPath.push(newPt);
    }
    newPath.push(path[path.length - 1]);

    return newPath;
  }

  /**
   * Translate the height map up or down so that a particular
   * quantile is at zero
   *
   * @param {object} landscape
   * @param {number} q
   * @returns {meshObject} newH
   */
  static setSeaLevel(landscape, q) {
    var newH = Mesh.zero(landscape.mesh),
      delta = Mesh.quantile(landscape, q);
    for (let i = 0; i < landscape.length; i++) {
      let landVal = landscape[i];
      newH[i] = landVal - delta;
    }

    return newH;
  }

  /**
   * Generate a constant slope think of it as tectonic
   * uplift on one side of the map
   *
   * @param {meshObject} mesh
   * @param {Array} direction
   */
  static slope(mesh, direction) {
    return mesh.map(function (x) {
      return x[0] * direction[0] + x[1] * direction[1];
    });
  }

  /**
   *
   * @param {object} landscape
   * @param {number} i
   * @returns {*[]}
   */
  static triSlope(landscape, i) {
    var nbs = Mesh.neighbors(landscape.mesh, i);
    if (nbs.length != 3) {
      return [0, 0];
    }
    var p0 = landscape.mesh.vxs[nbs[0]],
      p1 = landscape.mesh.vxs[nbs[1]],
      p2 = landscape.mesh.vxs[nbs[2]],

      x1 = p1[0] - p0[0],
      x2 = p2[0] - p0[0],
      y1 = p1[1] - p0[1],
      y2 = p2[1] - p0[1],

      det = x1 * y2 - x2 * y1,
      h1 = landscape[nbs[1]] - landscape[nbs[0]],
      h2 = landscape[nbs[2]] - landscape[nbs[0]];

    return [
      (y2 * h1 - y1 * h2) / det,
      (-x2 * h1 + x1 * h2) / det
    ];
  }

  /**
   *
   * @param {Array} pts
   * @param {extentObject} extent
   * @returns {d3.voronoi} v
   */
  static voronoi(pts, extent) {
    var w = extent.width / 2,
      h = extent.height / 2,
      v = d3.voronoi().extent([[-w, -h], [w, h]])(pts);

    return v;
  }

  /**
   *
   * @param {meshObject} mesh
   * @returns {meshObject} z
   */
  static zero(mesh) {
    var z = [];
    for (let i = 0; i < mesh.vxs.length; i++) {
      z[i] = 0;
    }
    z.mesh = mesh;

    return z;
  }

}

/**
 *
 */
class Map {

  /**
   *
   * @param {paramsObject} params
   * @param {ID3Selection} svg
   * @returns {Map}
   */
  constructor(params, svg) {
    this.settings = params || paramsObject;
    this.generator = this.settings.generator;
    this.svg = svg;

    this.borders = [];
    this.cities = [];
    this.coasts = [];
    this.rivers = [];
    this.territories = [];
    this.cityLabels = [];
    this.regionLabels = [];

    this.viewBorders = true;
    this.viewCities = true;
    this.viewCoasts = true;
    this.viewHeight = false;
    this.viewLabels = true;
    this.viewRivers = true;
    this.viewScores = false;
    this.viewSlopes = true;
    this.viewPoints = false;

    /**
     * Landscape functions
     */
    this.functions = {
      newEmptyMap: () => {
        this.clearMap();
        this.drawMap();
      },
      newCoastLine: () => {
        this.clearMap();
        this.landscape = generateCoastLine(this.settings);
        this.doMap();
      },
      newIsland: () => {
        this.clearMap();
        this.landscape = generateIsland(this.settings);
        this.doMap();
      },
      newLake: () => {
        this.clearMap();
        this.landscape = generateLake(this.settings);
        this.doMap();
      },
      /**
       * Landscape add actions
       */
      add: {
        city: () => {
          this.addCity();
          this.drawMap();
        },
        hill: () => {
          this.landscape = add(this.landscape, Mesh.cone(this.landscape.mesh, -0.5));
          this.drawMap();
        },
        lake: () => {
          this.landscape = add(this.landscape, Mesh.cone(this.landscape.mesh, 0.5));
          this.drawMap();
        },
        mountains: () => {
          this.landscape = add(this.landscape, Mesh.mountains(this.landscape.mesh, 5));
          this.drawMap();
        },
        randomSlope: () => {
          this.landscape = add(this.landscape, Mesh.slope(this.landscape.mesh, randomVector(4)));
          this.drawMap();
        }
      },
      /**
       * Landscape generate actions
       */
      generate: {
        borders: () => {
          this.getBorders(this.landscape, this.territories);
          this.drawMap();
        },
        cities: () => {
          this.getCities();
          this.drawMap();
        },
        coasts: () => {
          this.getCoasts(this.landscape);
          this.landscape = Mesh.setSeaLevel(this.landscape, this.settings.seaLevel);
          this.drawMap();
        },
        labels: () => {
          this.getLabels(this.landscape);
          this.drawLabels();
        },
        rivers: () => {
          this.getRivers(Mesh.zero(this.landscape.mesh), 0.01);
          this.drawMap();
        },
        territories: () => {
          this.getTerritories(this.landscape);
          this.drawMap();
        }
      },
      /**
       * Landscape actions
       */
      landscapeActions: {
        erode: () => {
          this.landscape = Mesh.doErosion(this.landscape);
          this.drawMap();
        },
        normalize: () => {
          this.landscape = Mesh.normalize(this.landscape);
          this.drawMap();
        },
        round: () => {
          this.landscape = Mesh.peaky(this.landscape);
          this.drawMap();
        },
        relax: () => {
          this.landscape = Mesh.relax(this.landscape);
          this.drawMap();
        },
        seaLevel: () => {
          this.landscape = Mesh.setSeaLevel(this.landscape, this.settings.seaLevel);
          this.drawMap();
        }
      },
      /**
       * Landscape toggle actions
       */
      toggle: {
        borders: () => {
          this.viewBorders = !this.viewBorders;
          this.drawMap();
        },
        cities: () => {
          this.viewCities = !this.viewCities;
          this.drawMap();
        },
        coasts: () => {
          this.viewCoasts = !this.viewCoasts;
          this.drawMap();
        },
        heightMap: () => {
          this.viewHeight = !this.viewHeight;
          this.drawMap();
        },
        labels: () => {
          this.viewLabels = !this.viewLabels;
          this.drawMap();
        },
        points: () => {
          this.viewPoints = !this.viewPoints;
          this.drawMap();
        },
        rivers: () => {
          this.viewRivers = !this.viewRivers;
          this.drawMap();
        },
        slopes: () => {
          this.viewSlopes = !this.viewSlopes;
          this.drawMap();
        }
      }
    };

    this.svg.attr('float', 'left').attr('background-color', 'white');

    return this;
  }

  /**
   *
   */
  addCity() {
    this.cities = this.cities || [];
    var score = this.cityScore(),
      newCity = d3.scan(score, d3.descending);
    this.cities.push(newCity);
  }

  /**
   * Generate a score for each point, which is a combination of three things:
   * Water flux - we want cities to be preferentially located on rivers, so high
   *  water flux gets a bonus
   * Distance from other cities - we want cities to be spread out, so penalize
   *  locations which are too close to an existing city
   * Distance from the edge of the map - the other two criteria alone tend to push
   *  cities to the map edge, which isn't ideal, so penalize locations too close to
   *  the edge
   *
   * @returns {object} score
   */
  cityScore() {
    var score = Mesh.map(Mesh.getFlux(this.landscape), Math.sqrt);
    for (let i = 0; i < this.landscape.length; i++) {
      if (this.landscape[i] <= 0 || Mesh.isNearEdge(this.landscape.mesh, i)) {
        score[i] = -999999;
        continue;
      }
      score[i] += 0.01 / (1e-9 + Math.abs(this.landscape.mesh.vxs[i][0]) - this.landscape.mesh.extent.width / 2);
      score[i] += 0.01 / (1e-9 + Math.abs(this.landscape.mesh.vxs[i][1]) - this.landscape.mesh.extent.height / 2);
      for (let j = 0; j < this.cities.length; j++) {
        score[i] -= 0.02 / (Mesh.distance(this.landscape.mesh, this.cities[j], i) + 1e-9);
      }
    }

    return score;
  }

  /**
   * Clear the map
   */
  clearMap() {
    this.borders = [];
    this.cities = [];
    this.coasts = [];
    this.rivers = [];
    this.territories = [];
    this.cityLabels = [];
    this.regionLabels = [];
    this.landscape = Mesh.zero(Mesh.generateGoodMesh(this.settings.numPts, this.settings.extent));

    return this;
  }

  /**
   * Generate all the coasts, borders, cities etc and draw them
   */
  doMap() {
    this.getRivers(this.landscape, 0.01);
    this.getCoasts(this.landscape);
    this.getCities();
    this.getTerritories(this.landscape);
    this.getBorders(this.landscape, this.territories);
    this.getLabels(this.landscape);

    this.drawMap();

    return this;
  }

  /**
   * Draw the cities
   *
   * @returns {Map}
   */
  drawCircles(cls, pts) {
    var circles = this.svg.selectAll('circle.' + cls).data(pts);
    circles.enter().append('circle').classed(cls, true);
    circles.exit().remove();

    this.svg.selectAll('circle.' + cls)
      .attr('cx', (d) => {
        return 1000 * this.landscape.mesh.vxs[d][0];
      })
      .attr('cy', (d) => {
        return 1000 * this.landscape.mesh.vxs[d][1];
      })
      .attr('r', (d, i) => {
        return i >= this.settings.numTerritories ? 4 : 10;
      })
      .style('fill', 'white')
      .style('stroke', 'orange')
      .style('stroke-width', 4)
      .style('stroke-linecap', 'round')
      .raise();

    return this;
  }

  /**
   * Draw the labels for regions and cities
   *
   * @param cls
   * @param pts
   */
  drawLabels(cls, pts) {
    var cityTexts = this.svg.selectAll('text.city').data(this.cityLabels);
    cityTexts.enter().append('text').classed('city', true);
    cityTexts.exit().remove();

    this.svg.selectAll('text')
      .style('font-family', '"Palatino Linotype", "Book Antiqua", "Palatino", "serif"')
      .style('color', 'black')
      .style('stroke', 'white')
      .style('stroke-linejoin', 'round')
      .style('paint-order', 'stroke');

    this.svg.selectAll('text.city')
      .attr('x', function (d) {
        return 1000 * d.x;
      })
      .attr('y', function (d) {
        return 1000 * d.y;
      })
      .style('font-size', function (d) {
        return d.size;
      })
      .style('text-anchor', function (d) {
        return d.align;
      })
      .text(function (d) {
        return d.text;
      })
      .style('stroke-width', 2)
      .raise();


    var regionTexts = this.svg.selectAll('text.region').data(this.regionLabels);
    regionTexts.enter().append('text').classed('region', true);
    regionTexts.exit().remove();

    this.svg.selectAll('text.region')
      .attr('x', function (d) {
        return 1000 * d.x;
      })
      .attr('y', function (d) {
        return 1000 * d.y;
      })
      .style('font-size', function (d) {
        return 1000 * d.size;
      })
      .text(function (d) {
        return d.text;
      })
      .style('font-variant', 'small-caps')
      .style('text-anchor', 'middle')
      .style('stroke-width', 5)
      .raise();

    return this;
  }

  /**
   * Draw the map and all the details
   *
   * @returns {Map}
   */
  drawMap() {
    this.drawPaths('river', this.viewRivers ? this.rivers : []);
    this.drawPaths('coast', this.viewCoasts ? this.coasts : []);
    this.drawPaths('border', this.viewBorders ? this.borders : []);
    this.drawSlopes(this.viewSlopes ? this.landscape : Mesh.zero(this.landscape.mesh));
    this.drawCircles('city', this.viewCities ? this.cities : []);
    this.drawLabels(this.viewLabels ? this.cityLabels : []);
    this.drawLabels(this.viewLabels ? this.regionLabels : []);

    if (this.viewPoints) {
      this.drawPoints(this.viewPoints ? this.landscape.mesh.pts : []);
    }
    if (this.viewHeight) {
      this.drawVoronoi(this.landscape, -1, 1);
    }

    if (this.viewErosionRate) {
      this.drawVoronoi(Mesh.erosionRate(this.landscape));
    }

    if (this.viewScores) {
      this.drawScores();
    }

    return this;
  }

  /**
   * Draw the path of points for the passed in class, river, coast etc.
   *
   * @param {string} cls
   * @param {Array} paths
   * @returns {Map}
   */
  drawPaths(cls, pts) {
    var svgPaths = this.svg.selectAll('path.' + cls).data(pts);
    svgPaths.enter().append('path').classed(cls, true);
    svgPaths.exit().remove();

    this.svg.selectAll('path.' + cls)
      .attr('d', this.makeD3Path);

    this.svg.selectAll('path, line')
      .style('fill', 'none')
      .style('stroke', 'black')
      .style('stroke-linecap', 'round');

    this.svg.selectAll('path.field')
      .style('stroke', 'none')
      .style('fill-opacity', 1.0);

    this.svg.selectAll('path.border')
      .style('stroke', 'red')
      .style('stroke-width', 5)
      .style('stroke-dasharray', '4,4')
      .style('stroke-linecap', 'butt');

    this.svg.selectAll('path.coast')
      .style('stroke-width', 3);

    this.svg.selectAll('path.river')
      .style('stroke', 'blue')
      .style('stroke-width', 2);

    return this;
  }

  /**
   * Draw the points of the map, mainly a way to see the underlying mesh points
   *
   * @returns {Map}
   */
  drawPoints(pts) {
    var circle = this.svg.selectAll('circle').data(pts);
    circle.enter().append('circle');
    circle.exit().remove();
    d3.selectAll('circle')
      .attr('cx', function (d) {
        return 1000 * d[0];
      })
      .attr('cy', function (d) {
        return 1000 * d[1];
      })
      .attr('r', 100 / Math.sqrt(pts.length));

    return this;
  }

  /**
   * Draw teh scores of the cities
   */
  drawScores() {
    var score = this.cityScore();
    this.drawVoronoi(score, d3.max(score) - 0.5, d3.max(score) + 0.5);

    return this;
  }

  /**
   * Draw strokes which go up and right if the terrain slopes upwards from
   * left to right, and down and right if the terrain slopes downwards.
   * Similarly, the strokes on the 'near' side of hills should be longer than
   * those on the 'far' side.
   *
   * @returns {Map}
   */
  drawSlopes(landscape) {
    var strokes = [],
      r = 0.25 / Math.sqrt(landscape.length);
    for (let i = 0; i < landscape.length; i++) {
      if (landscape[i] <= 0 || Mesh.isNearEdge(landscape.mesh, i)) {
        continue;
      }
      var nbs = Mesh.neighbors(landscape.mesh, i),
        s = 0,
        s2 = 0;
      nbs.push(i);
      for (let j = 0; j < nbs.length; j++) {
        var slopes = Mesh.triSlope(landscape, nbs[j]);
        s += slopes[0] / 10;
        s2 += slopes[1];
      }
      s /= nbs.length;
      s2 /= nbs.length;
      if (Math.abs(s) < runIf(0.1, 0.4)) {
        continue;
      }
      var l = r * runIf(1, 2) * (1 - 0.2 * Math.pow(Math.atan(s), 2)) * Math.exp(s2 / 100),
        x = landscape.mesh.vxs[i][0],
        y = landscape.mesh.vxs[i][1];
      if (Math.abs(l * s) > 2 * r) {
        var n = Math.floor(Math.abs(l * s / r));
        l /= n;
        if (n > 4) {
          n = 4;
        }
        for (let j = 0; j < n; j++) {
          var u = rNorm() * r,
            v = rNorm() * r;
          strokes.push([[x + u - l, y + v + l * s], [x + u + l, y + v - l * s]]);
        }
      } else {
        strokes.push([[x - l, y + l * s], [x + l, y - l * s]]);
      }
    }

    var lines = this.svg.selectAll('line.slope').data(strokes);
    lines.enter().append('line').classed('slope', true);
    lines.exit().remove();

    this.svg.selectAll('line.slope')
      .attr('x1', function (d) {
        return 1000 * d[0][0];
      })
      .attr('y1', function (d) {
        return 1000 * d[0][1];
      })
      .attr('x2', function (d) {
        return 1000 * d[1][0];
      })
      .attr('y2', function (d) {
        return 1000 * d[1][1];
      });

    this.svg.selectAll('line.slope')
      .style('stroke', 'green')
      .style('stroke-width', 1);

    return this;
  }

  /**
   * Draw the points as triangles and colors indicating height
   *
   * @param {Object} landscape
   * @param {number|null} lo
   * @param {number|null} hi
   */
  drawVoronoi(landscape, lo = null, hi = null) {
    if (hi == null) {
      hi = d3.max(landscape) + 1e-9;
    }
    if (lo == null) {
      lo = d3.min(landscape) - 1e-9;
    }
    var mappedVals = landscape.map(function (x) {
        return x > hi ? 1 : x < lo ? 0 : (x - lo) / (hi - lo);
      }),
      tris = this.svg.selectAll('path.field').data(landscape.mesh.tris);

    tris.enter().append('path').classed('field', true);
    tris.exit().remove();

    this.svg.selectAll('path.field')
      .attr('d', this.makeD3Path)
      .style('fill', function (d, i) {
        return d3.interpolateViridis(mappedVals[i]);
      });

    return this;
  }

  /**
   * Generate and set the borders for the landscape
   *
   * @returns {Map}
   */
  getBorders(landscape, territories) {
    var edges = [];
    for (let i = 0; i < territories.mesh.edges.length; i++) {
      var e = territories.mesh.edges[i];
      if (e[3] == undefined) {
        continue;
      }
      if (Mesh.isNearEdge(territories.mesh, e[0]) || Mesh.isNearEdge(territories.mesh, e[1])) {
        continue;
      }
      if (landscape[e[0]] < 0 || landscape[e[1]] < 0) {
        continue;
      }
      if (territories[e[0]] != territories[e[1]]) {
        edges.push([e[2], e[3]]);
      }
    }

    this.borders = Mesh.mergeSegments(edges).map(Mesh.relaxPath);

    return this;
  }

  /**
   * Generate and set the cities for the landscape
   *
   * @returns {Map}
   */
  getCities() {
    var n = this.settings.numCities;
    for (let i = 0; i < n; i++) {
      var score = this.cityScore(),
        newCity = d3.scan(score, d3.descending);
      this.cities.push(newCity);
    }

    return this;
  }

  /**
   * Generate and set the coasts for the landscape
   *
   * @returns {Map}
   */
  getCoasts(landscape) {
    this.coasts = Mesh.contour(landscape, 0);

    return this;
  }

  /**
   * Generate and set the labels for the landscape
   *
   * @returns {Map}
   */
  getLabels(landscape) {
    var numTer = this.settings.numTerritories,
      avoids = [this.rivers, this.coasts, this.borders],
      lang = makeRandomLanguage();
    this.cityLabels = [];

    var penalty = (label) => {
      let pen = 0;
      if (label.x0 < -0.45 * landscape.mesh.extent.width) {
        pen += 100;
      }
      if (label.x1 > 0.45 * landscape.mesh.extent.width) {
        pen += 100;
      }
      if (label.y0 < -0.45 * landscape.mesh.extent.height) {
        pen += 100;
      }
      if (label.y1 > 0.45 * landscape.mesh.extent.height) {
        pen += 100;
      }

      for (let i = 0; i < this.cityLabels.length; i++) {
        let oLabel = this.cityLabels[i];
        if (label.x0 < oLabel.x1 && label.x1 > oLabel.x0 &&
          label.y0 < oLabel.y1 && label.y1 > oLabel.y0) {
          pen += 100;
        }
      }

      for (let i = 0; i < this.cities.length; i++) {
        let c = landscape.mesh.vxs[this.cities[i]];
        if (label.x0 < c[0] && label.x1 > c[0] && label.y0 < c[1] && label.y1 > c[1]) {
          pen += 100;
        }
      }

      for (let i = 0; i < avoids.length; i++) {
        let avoid = avoids[i];
        for (let j = 0; j < avoid.length; j++) {
          let avpath = avoid[j];
          for (let k = 0; k < avpath.length; k++) {
            let pt = avpath[k];
            if (pt[0] > label.x0 && pt[0] < label.x1 && pt[1] > label.y0 && pt[1] < label.y1) {
              pen++;
            }
          }
        }
      }

      return pen;
    };

    // City labels
    for (let i = 0; i < this.cities.length; i++) {
      let x = landscape.mesh.vxs[this.cities[i]][0],
        y = landscape.mesh.vxs[this.cities[i]][1],
        text = makeName(lang, 'city'),
        town = i < numTer ? false : true,
        size = town ? this.settings.fontSizes.city : this.settings.fontSizes.town,
        sx = 0.65 * size / 1000 * text.length,
        sy = size / 1000,
        possLabels = [
          {
            x: x + 0.8 * sy,
            y: y + 0.3 * sy,
            align: 'start',
            x0: x + 0.7 * sy,
            y0: y - 0.6 * sy,
            x1: x + 0.7 * sy + sx,
            y1: y + 0.6 * sy
          },
          {
            x: x - 0.8 * sy,
            y: y + 0.3 * sy,
            align: 'end',
            x0: x - 0.9 * sy - sx,
            y0: y - 0.7 * sy,
            x1: x - 0.9 * sy,
            y1: y + 0.7 * sy
          },
          {
            x: x,
            y: y - 0.8 * sy,
            align: 'middle',
            x0: x - sx / 2,
            y0: y - 1.9 * sy,
            x1: x + sx / 2,
            y1: y - 0.7 * sy
          },
          {
            x: x,
            y: y + 1.2 * sy,
            align: 'middle',
            x0: x - sx / 2,
            y0: y + 0.1 * sy,
            x1: x + sx / 2,
            y1: y + 1.3 * sy
          }
        ],
        label = possLabels[d3.scan(possLabels, function (a, b) {
          return penalty(a) - penalty(b);
        })];
      label.text = text;
      label.size = size;
      this.cityLabels.push(label);
    }

    // Region labels
    this.regionLabels = [];
    for (let i = 0; i < numTer; i++) {
      let city = this.cities[i],
        text = makeName(lang, 'region'),
        sy = this.settings.fontSizes.region / 1000,
        sx = 0.6 * text.length * sy,
        lc = this.getTerritoryCenter(this.territories, city, true),
        oc = this.getTerritoryCenter(this.territories, city, false),
        best = 0,
        bestScore = -999999;
      for (let j = 0; j < landscape.length; j++) {
        let score = 0,
          v = landscape.mesh.vxs[j];
        score -= 3000 * Math.sqrt((v[0] - lc[0]) * (v[0] - lc[0]) + (v[1] - lc[1]) * (v[1] - lc[1]));
        score -= 1000 * Math.sqrt((v[0] - oc[0]) * (v[0] - oc[0]) + (v[1] - oc[1]) * (v[1] - oc[1]));
        if (this.territories[j] != city) {
          score -= 3000;
        }
        for (let k = 0; k < this.cities.length; k++) {
          let u = landscape.mesh.vxs[this.cities[k]];
          if (Math.abs(v[0] - u[0]) < sx && Math.abs(v[1] - sy / 2 - u[1]) < sy) {
            score -= k < numTer ? 4000 : 500;
          }
          if (v[0] - sx / 2 < this.cityLabels[k].x1 &&
            v[0] + sx / 2 > this.cityLabels[k].x0 &&
            v[1] - sy < this.cityLabels[k].y1 &&
            v[1] > this.cityLabels[k].y0) {
            score -= 5000;
          }
        }
        for (let k = 0; k < this.regionLabels.length; k++) {
          let label = this.regionLabels[k];
          if (v[0] - sx / 2 < label.x + label.width / 2 &&
            v[0] + sx / 2 > label.x - label.width / 2 &&
            v[1] - sy < label.y &&
            v[1] > label.y - label.size) {
            score -= 20000;
          }
        }
        if (landscape[j] <= 0) {
          score -= 500;
        }
        if (v[0] + sx / 2 > 0.5 * landscape.mesh.extent.width) {
          score -= 50000;
        }
        if (v[0] - sx / 2 < -0.5 * landscape.mesh.extent.width) {
          score -= 50000;
        }
        if (v[1] > 0.5 * landscape.mesh.extent.height) {
          score -= 50000;
        }
        if (v[1] - sy < -0.5 * landscape.mesh.extent.height) {
          score -= 50000;
        }
        if (score > bestScore) {
          bestScore = score;
          best = j;
        }
      }
      this.regionLabels.push({
        text: text,
        x: landscape.mesh.vxs[best][0],
        y: landscape.mesh.vxs[best][1],
        size: sy,
        width: sx
      });
    }

    return this;
  }

  /**
   * Generate and set the rivers for the landscape
   *
   * @param {number} limit
   * @returns {Map}
   */
  getRivers(landscape, limit = 0.01) {
    var dh = Mesh.downhill(landscape),
      flux = Mesh.getFlux(landscape),
      links = [],
      above = 0;
    this.rivers = [];

    for (let i = 0; i < landscape.length; i++) {
      if (landscape[i] > 0) {
        above++;
      }
    }

    limit *= above / landscape.length;
    for (let i = 0; i < dh.length; i++) {
      if (Mesh.isNearEdge(landscape.mesh, i)) {
        continue;
      }
      if (flux[i] > limit && landscape[i] > 0 && dh[i] >= 0) {
        let up = landscape.mesh.vxs[i],
          down = landscape.mesh.vxs[dh[i]];
        if (landscape[dh[i]] > 0) {
          links.push([up, down]);
        } else {
          links.push([up, [(up[0] + down[0]) / 2, (up[1] + down[1]) / 2]]);
        }
      }
    }

    this.rivers = Mesh.mergeSegments(links).map(Mesh.relaxPath);

    return this;
  }

  /**
   * Expand regions outwards from each city, so that each region consists of the
   * points which are 'closest' to its city, according to a particular distance
   * measure. This distance measure is calculated by adding up the cost of the route,
   * based on these criteria:
   * Horizontal distance
   * Slope - uphill is much cheaper than downhill, so regions expand until they hit
   *  the top of ridges, then stop
   * Water flux - crossing a river is expensive
   * Shorelines - there is a large penalty for going from land to water (or vice versa)
   *  and a smaller penalty for travelling by water
   *
   * @returns {Map}
   */
  getTerritories(landscape) {
    this.territories = [];
    var n = this.settings.numTerritories;
    if (n > this.cities.length) {
      n = this.cities.length;
    }
    let flux = Mesh.getFlux(landscape),
      queue = new PriorityQueue({
        comparator: function (a, b) {
          return a.score - b.score;
        }
      });

    var weight = (u, v) => {
      var horiz = Mesh.distance(landscape.mesh, u, v),
        vert = landscape[v] - landscape[u];
      if (vert > 0) {
        vert /= 10;
      }
      let diff = 1 + 0.25 * Math.pow(vert / horiz, 2);
      diff += 100 * Math.sqrt(flux[u]);
      if (landscape[u] <= 0) {
        diff = 100;
      }
      if ((landscape[u] > 0) != (landscape[v] > 0)) {
        return 1000;
      }

      return horiz * diff;
    };

    for (let i = 0; i < n; i++) {
      this.territories[this.cities[i]] = this.cities[i];
      let nbs = Mesh.neighbors(landscape.mesh, this.cities[i]);
      for (let j = 0; j < nbs.length; j++) {
        queue.queue({
          score: weight(this.cities[i], nbs[j]),
          city: this.cities[i],
          vx: nbs[j]
        });
      }
    }

    while (queue.length) {
      let u = queue.dequeue();
      if (this.territories[u.vx] != undefined) {
        continue;
      }
      this.territories[u.vx] = u.city;
      var nbs = Mesh.neighbors(landscape.mesh, u.vx);
      for (let i = 0; i < nbs.length; i++) {
        let v = nbs[i],
          newDist = weight(u.vx, v);

        if (this.territories[v] != undefined) {
          continue;
        }
        queue.queue({
          score: u.score + newDist,
          city: u.city,
          vx: v
        });
      }
    }
    this.territories.mesh = landscape.mesh;

    return this;
  }

  /**
   * Find the center of a territory
   *
   * @param {Object} terr
   * @param {Object} city
   * @param {boolean} landOnly
   * @returns {*[]}
   */
  getTerritoryCenter(terr, city, landOnly) {
    var x = 0, y = 0, n = 0;
    for (let i = 0; i < terr.length; i++) {
      if (terr[i] != city) {
        continue;
      }
      if (landOnly && this.landscape[i] <= 0) {
        continue;
      }
      x += terr.mesh.vxs[i][0];
      y += terr.mesh.vxs[i][1];
      n++;
    }

    return [x / n, y / n];
  }

  /**
   * Make a set of d3 coordinates for a path
   *
   * @param {Array} path
   * @return {string}
   */
  makeD3Path(path) {
    var p = d3.path();
    p.moveTo(1000 * path[0][0], 1000 * path[0][1]);
    for (let i = 1; i < path.length; i++) {
      p.lineTo(1000 * path[i][0], 1000 * path[i][1]);
    }

    return p.toString();
  }

}
