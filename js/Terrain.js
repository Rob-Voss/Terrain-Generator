"use strict";

let perlin = noise,

  /**
   *
   */
  rNorm = (function () {
    let z2 = null;

    function rNorm() {
      if (z2 != null) {
        let tmp = z2;
        z2 = null;
        return tmp;
      }
      let x1 = 0,
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

    return rNorm;
  })();

/**
 * Default font sizes for the map
 * @typedef {object} fontSizeObject
 * @property {number} region
 * @property {number} city
 * @property {number} town
 */
const fontSizeObject = {
    region: 40,
    city: 25,
    town: 20
  },
  /**
   * Extent Object
   * @typedef {object} extentObject
   * @property {number} width
   * @property {number} height
   */
  extentObject = {
    width: 1,
    height: 1
  },
  /**
   * Scale Object
   * @typedef {object} scaleObject
   * @property {number} width
   * @property {number} height
   */
  scaleObject = {
    width: 1000,
    height: 1000
  },
  /**
   * Default options for the map generator
   * @typedef {object} paramsObject
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
   * @typeDef {object} meshObject
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
  },
  /**
   * Landscape Object
   * @typeDef {Array} landscapeObject
   * @property {meshObject} mesh
   */
  landscapeObject = [{
    mesh: meshObject
  }];

/**
 *
 * @param {number} gradIntensity
 * @param {number} perlIntensity
 */
function generateHex(svg, gradIntensity = 1, perlIntensity = 50) {
  let l = 100,
    grid = Hex.getGrid(l),
    hl = Math.floor(l / 2),
    center = grid[hl][hl],
    circles = [],
    color = 0,
    width = 10,
    height = width * Math.sqrt(3) / 2,
    radius = width / 2;

  for (let i = 1; i <= hl; i += 1) {
    circles[i] = Hex.getRing(center, i);
  }

  let gradient = Hex.getGradientColors(grid, circles, gradIntensity),
    perlinC = Hex.getPerlinColors(grid, perlIntensity);

  for (let i = 0; i < grid.length; i += 1) {
    for (let j = 0; j < grid[i].length; j += 1) {
      let x = width * 3 / 4 * (i + 1),
        y = (i % 2 !== 0) ? (radius * Math.sqrt(3) / 2 + height * j) : height * j,
        colorValue = perlinC[i][j] + gradient[i][j];

      if (colorValue < 50) {
        color = '#567A32';
      } else if (colorValue < 100) {
        color = '#78ab46';
      } else if (colorValue < 150) {
        color = '#bced91';
      } else if (colorValue < 180) {
        color = '#dbd1b4';
      } else {
        color = '#b4d2db';
      }

      svg.select('g')
        .append('polygon')
        .attr('points', Hex.getHexagonPoints(x, y, radius))
        .attr('fill', color)
        .attr('stroke', 'grey')
        .attr('stroke-width', 0.5);
    }
  }
}

/**
 *
 * @param {paramsObject} params
 * @return {landscapeObject} landscape
 */
function generateCoastLine(params = paramsObject) {
  let mesh = Mesh.generateGoodMesh(params.numPts, params.extent),
    landscape = add(
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
 * @param {paramsObject} params
 */
function generateLake(params = paramsObject) {
  let mesh = Mesh.generateGoodMesh(params.numPts, params.extent),

    landscape = add(
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
 * @param {paramsObject} params
 */
function generateIslands(params = paramsObject) {
  let mesh = Mesh.generateGoodMesh(params.numPts, params.extent),

    landscape = add(
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

/**
 *
 * @param {number} scale
 * @return {[*,*]}
 */
function randomVector(scale) {
  return [scale * rNorm(), scale * rNorm()];
}

/**
 *
 * @param {number} lo
 * @param {number} hi
 * @return {*}
 */
function runIf(lo, hi) {
  return lo + Math.random() * (hi - lo);
}

/**
 *
 * @return {landscapeObject}
 */
function add() {
  let n = arguments[0].length,
    newVals = Mesh.zero(arguments[0].mesh);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < arguments.length; j++) {
      newVals[i] += arguments[j][i];
    }
  }

  return newVals;
}

/**
 * @class Hex
 */
class Hex {

  /**
   *
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {Hex}
   */
  constructor(x = 0, y = 0, z = 0) {
    this.x = x;
    this.y = y;
    this.z = z;

    return this;
  }

  /**
   *
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {Hex}
   */
  static create(x, y, z) {
    return new Hex(x, y, z);
  }

  /**
   *
   * @param {Hex} a
   * @param {Hex} b
   * @return {*|Hex}
   */
  static add(a, b) {
    return Hex.create(a.x + b.x, a.y + b.y, a.z + b.z);
  }

  /**
   *
   * @param {Hex} a
   * @param {number} k
   * @return {*|Hex}
   */
  static scale(a, k) {
    return Hex.create(a.x * k, a.y * k, a.z * k);
  }

  /**
   *
   * @param {Hex} hex
   * @param {number} direction
   * @return {*|Hex}
   */
  static getNeighbour(hex, direction) {
    return Hex.add(hex, Hex.getDirection(direction));
  }

  /**
   *
   * @param {number} direction
   * @return {*}
   */
  static getDirection(direction) {
    let directions = [
      Hex.create(1, 0, -1),
      Hex.create(1, -1, 0),
      Hex.create(0, -1, 1),
      Hex.create(-1, 0, 1),
      Hex.create(-1, 1, 0),
      Hex.create(0, 1, -1)
    ];

    return directions[direction];
  }

  /**
   *
   * @param {Hex} hex
   * @return {[Hex]}
   */
  static getAllNeighbours(hex) {
    return [
      Hex.getNeighbour(hex, 0),
      Hex.getNeighbour(hex, 1),
      Hex.getNeighbour(hex, 2),
      Hex.getNeighbour(hex, 3),
      Hex.getNeighbour(hex, 4),
      Hex.getNeighbour(hex, 5)
    ];
  }

  /**
   *
   * @param {number} col
   * @param {number} row
   * @return {*}
   */
  static convertOffsetToCube(col, row) {
    var x = col,
      z = (col % 2 === 0) ? row - (col + 1) / 2 : row - col / 2,
      y = -x - z;

    return Hex.round(Hex.create(x, y, z));
  }

  /**
   *
   * @param {Hex} hex
   * @return {*|Hex}
   */
  static round(hex) {
    var x = Math.trunc(Math.round(hex.x)),
      y = Math.trunc(Math.round(hex.y)),
      z = Math.trunc(Math.round(hex.z)),
      xDiff = Math.abs(x - hex.x),
      yDiff = Math.abs(y - hex.y),
      zDiff = Math.abs(z - hex.z);

    if (xDiff > yDiff && xDiff > zDiff) {
      x = -y - z;
    } else {
      if (yDiff > zDiff) {
        y = -x - z;
      } else {
        z = -x - y;
      }
    }
    return Hex.create(x, y, z);
  }

  /**
   *
   * @param {Hex} a
   * @param {Hex} b
   * @return {Array}
   */
  static getLinedraw(a, b) {
    let distance = Hex.getDistance(a, b),
      results = [],
      step = 1.0 / Math.max(distance, 1),
      i;

    for (i = 0; i <= distance; i += 1) {
      results.push(Hex.round(Hex.lerp(a, b, step * i)));
    }

    return results;
  }

  /**
   *
   * @param {Hex} a
   * @param {Hex} b
   * @return {*}
   */
  static getDistance(a, b) {
    return Hex.getLength(Hex.subtract(a, b));
  }

  /**
   *
   * @param {Hex} a
   * @param {Hex} b
   * @return {*|Hex}
   */
  static subtract(a, b) {
    return Hex.create(a.x - b.x, a.y - b.y, a.z - b.z);
  }

  /**
   *
   * @param {Hex} hex
   * @return {number}
   */
  static getLength(hex) {
    return Math.trunc((Math.abs(hex.x) + Math.abs(hex.y) + Math.abs(hex.z)) / 2);
  }

  /**
   *
   * @param {Hex} a
   * @param {Hex} b
   * @param {number} t
   * @return {*|Hex}
   */
  static lerp(a, b, t) {
    return Hex.create(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, a.z + (b.z - a.z) * t);
  }

  /**
   *
   * @param {number} center
   * @param {number} radius
   * @return {Array}
   */
  static getRing(center, radius) {
    let results = [],
      cube = Hex.add(center, Hex.scale(Hex.create(-1, 1, 0), radius));
    for (let i = 0; i < 6; i += 1) {
      for (let j = 0; j < radius; j += 1) {
        results.push(cube);
        cube = Hex.getNeighbour(cube, i);
      }
    }
    return results;
  }

  /**
   *
   * @param {number} l
   * @return {Array}
   */
  static getGrid(l) {
    let grid = [];
    for (let i = 0; i < l; i += 1) {
      grid[i] = [];
      for (let j = 0; j < l; j += 1) {
        grid[i][j] = Hex.convertOffsetToCube(i, j);
      }
    }
    return grid;
  }

  /**
   *
   * @param {Array} grid
   * @param {Array} circles
   * @param {number} intensity
   * @return {Array}
   */
  static getGradientColors(grid, circles, intensity = 1) {
    let colors = [];

    for (let i = 0; i < grid.length; i += 1) {
      colors[i] = [];
      for (let j = 0; j < grid[i].length; j += 1) {
        let color = 255 * intensity,
          iZ = i === Math.floor(grid.length / 2),
          jZ = j === Math.floor(grid[i].length / 2);

        if (iZ && jZ) {
          color = 0;
        } else {
          for (let ci = 1; ci < circles.length; ci += 1) {
            for (let cj = 0; cj < circles[ci].length; cj += 1) {
              let cir = circles[ci][cj],
                grd = grid[i][j];
              if (cir.x === grd.x && cir.y === grd.y && cir.z === grd.z) {
                color = Math.round(255 - (circles.length - ci) * 255 / circles.length) * intensity;
              }
            }
          }
        }
        colors[i][j] = color;
      }
    }
    return colors;
  }

  /**
   *
   * @param {Array} grid
   * @param {number} intensity
   * @return {Array}
   */
  static getPerlinColors(grid, intensity = 50) {
    let colors = [];
    noise.seed(Math.random());

    for (let i = 0; i < grid.length; i += 1) {
      colors[i] = [];
      for (let j = 0; j < grid[i].length; j += 1) {
        let value = noise.simplex2(i / intensity, j / intensity);
        value *= 256;

        colors[i][j] = Math.round(Math.abs(value));
      }
    }

    return colors;
  }

  /**
   *
   * @param {number} x
   * @param {number} y
   * @param {number} radius
   * @return {Array}
   */
  static getHexagonPoints(x, y, radius) {
    let points = [];
    for (let i = 0; i < 6; i += 1) {
      let angleDeg = 60 * i,
        angleRad = Math.PI / 180 * angleDeg,
        pointX = x + radius * Math.cos(angleRad),
        pointY = y + radius * Math.sin(angleRad);
      points.push([pointX, pointY]);
    }

    return points;
  }
};

/**
 * @class Mesh
 */
class Mesh {

  /**
   * The point where the three medians of the triangle intersect
   * @param {Array} pts
   * @return {*[]}
   */
  static centroid(pts) {
    let x = 0, y = 0;
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
   * @param {landscapeObject} landscape
   * @param {number} iterations
   * @return {*}
   */
  static cleanCoast(landscape, iterations) {
    for (let it = 0; it < iterations; it++) {
      let changed = 0,
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
   * @param {landscapeObject} landscape
   * @param {number} level
   * @return {Array}
   */
  static contour(landscape, level = 0) {
    let edges = [];
    for (let i = 0; i < landscape.mesh.edges.length; i++) {
      let e = landscape.mesh.edges[i],
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
   * @param {meshObject} mesh
   * @param {number} i
   * @param {number} j
   * @return {number}
   */
  static distance(mesh, i, j) {
    let p = mesh.vxs[i],
      q = mesh.vxs[j];

    return Math.sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]));
  }

  /**
   * Erode the coast line of the landscape
   * @param {landscapeObject} landscape
   * @param {number} amount
   * @param {number} n
   * @return {landscapeObject} landscape
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
   * @param {landscapeObject} landscape
   * @return {Array} downs
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

    let downs = [];
    for (let i = 0; i < landscape.length; i++) {
      downs[i] = downFrom(i);
    }
    landscape.downhill = downs;

    return downs;
  }

  /**
   *
   * @param {landscapeObject} landscape
   * @param {number} p
   * @return {meshObject} newH
   */
  static dropEdge(landscape, p) {
    p = p || 4;
    let newH = Mesh.zero(landscape.mesh);
    for (let i = 0; i < landscape.length; i++) {
      let v = landscape.mesh.vxs[i],
        x = 2.4 * v[0] / landscape.mesh.extent.width,
        y = 2.4 * v[1] / landscape.mesh.extent.height;
      newH[i] = landscape[i] - Math.exp(10 * (Math.pow(Math.pow(x, p) + Math.pow(y, p), 1 / p) - 1));
    }

    return newH;
  }

  /**
   * Generate a new landscape with eroded coast lines?
   * @param {landscapeObject} landscape
   * @param {number} amount
   * @return {meshObject} newH
   */
  static erode(landscape, amount) {
    let er = Mesh.erosionRate(landscape),
      newH = Mesh.zero(landscape.mesh),
      maxR = d3.max(er);
    for (let i = 0; i < landscape.length; i++) {
      newH[i] = landscape[i] - amount * (er[i] / maxR);
    }

    return newH;
  }

  /**
   * Set and erode the landscape rate?
   * @param {landscapeObject} landscape
   * @return {meshObject} newH
   */
  static erosionRate(landscape) {
    let flux = Mesh.getFlux(landscape),
      slope = Mesh.getSlope(landscape),
      newH = Mesh.zero(landscape.mesh);
    for (let i = 0; i < landscape.length; i++) {
      let river = Math.sqrt(flux[i]) * slope[i],
        creep = slope[i] * slope[i],
        total = 1000 * river + creep;
      total = total > 200 ? 200 : total;
      newH[i] = total;
    }

    return newH;
  }

  /**
   * Planchon-Darboux algorithm
   * @param {landscapeObject} landscape
   * @param {number} epsilon
   * @return {meshObject} newH
   */
  static fillSinks(landscape, epsilon = 1e-5) {
    let infinity = 999999,
      newH = Mesh.zero(landscape.mesh);
    for (let i = 0; i < landscape.length; i++) {
      if (Mesh.isNearEdge(landscape.mesh, i)) {
        newH[i] = landscape[i];
      } else {
        newH[i] = infinity;
      }
    }

    while (true) {
      let changed = false;
      for (let i = 0; i < landscape.length; i++) {
        if (newH[i] == landscape[i]) {
          continue;
        }
        let nbs = Mesh.neighbors(landscape.mesh, i);
        for (let j = 0; j < nbs.length; j++) {
          if (landscape[i] >= newH[nbs[j]] + epsilon) {
            newH[i] = landscape[i];
            changed = true;
            break;
          }
          let oh = newH[nbs[j]] + epsilon;
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
   * @param {landscapeObject} landscape
   */
  static findSinks(landscape) {
    let dh = Mesh.downhill(landscape),
      sinks = [];
    for (let i = 0; i < dh.length; i++) {
      let node = i;
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
   * @return {meshObject}
   */
  static generateGoodMesh(n = 1, extent = extentObject) {
    let pts = Mesh.generateGoodPoints(n, extent);

    return Mesh.makeMesh(pts, extent);
  }

  /**
   *
   * @param {number} n
   * @param {extentObject} extent
   * @return {Array}
   */
  static generateGoodPoints(n = 1, extent = extentObject) {
    // let pts = Mesh.generatePerlinNoise(n, extent, {width: 1000, height: 1000});
    let pts = Mesh.generatePoints(n, extent);
    pts = pts.sort(function (a, b) {
      return a[0] - b[0];
    });

    return Mesh.improvePoints(pts, 1, extent);
  }

  /**
   *
   * @param n
   * @param extent
   * @param scale
   * @return {Array}
   */
  static generatePerlinNoise(n = 1, extent = extentObject, scale = scaleObject) {
    let pts = [],
      max = -Infinity,
      min = Infinity,
      w = scale.width / 2,
      h = scale.height / 2;
    for (let x = 0; x < w; x++) {
      for (let y = 0; y < h; y++) {
        let value = perlin.simplex3(w / 100, h / 100, extent.height);
        if (max < value) {
          max = value;
        }
        if (min > value) {
          min = value;
        }
        let xP = (Math.random() - value) * extent.width,
          yP = (Math.random() - value) * extent.height;
        pts.push([xP, yP]);
      }
    }
    return pts;
  }

  /**
   *
   * @param {number} n
   * @param {extentObject} extent
   * @return {Array} pts
   */
  static generatePoints(n = 1, extent = extentObject) {
    let pts = [];
    for (let i = 0; i < n; i++) {
      let x = (Math.random() - 0.5),
        y = (Math.random() - 0.5);
      pts.push([
        x * extent.width,
        y * extent.height
      ]);
    }

    return pts;
  }

  /**
   * Get the water flux from each point to its 'downhill point'.
   * @param {landscapeObject} landscape
   * @return {meshObject} flux
   */
  static getFlux(landscape) {
    let dh = Mesh.downhill(landscape),
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
      let j = idxs[i];
      if (dh[j] >= 0) {
        flux[dh[j]] += flux[j];
      }
    }

    return flux;
  }

  /**
   * Get the slope values
   * @param {landscapeObject} landscape
   * @return {meshObject} slope
   */
  static getSlope(landscape) {
    let dh = Mesh.downhill(landscape),
      slope = Mesh.zero(landscape.mesh);
    for (let i = 0; i < landscape.length; i++) {
      let s = Mesh.triSlope(landscape, i);
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
   * @return {Array} pts
   */
  static improvePoints(pts, n = 1, extent = extentObject) {
    for (let i = 0; i < n; i++) {
      pts = Mesh.voronoi(pts, extent).polygons(pts).map(Mesh.centroid);
    }

    return pts;
  }

  /**
   * Check if it is an edge
   * @param {meshObject} mesh
   * @param {number} i
   * @return {boolean}
   */
  static isEdge(mesh, i) {
    return (mesh.adj[i].length < 3);
  }

  /**
   * Check if it is near an edge
   * @param {meshObject} mesh
   * @param {number} i
   * @return {boolean}
   */
  static isNearEdge(mesh, i) {
    let x = mesh.vxs[i][0],
      y = mesh.vxs[i][1],
      w = mesh.extent.width,
      h = mesh.extent.height;

    return x < -0.45 * w || x > 0.45 * w || y < -0.45 * h || y > 0.45 * h;
  }

  /**
   *
   * @param {Array} points
   * @param {extentObject} extent
   * @return {meshObject} mesh
   */
  static makeMesh(points, extent) {
    let voronoi = Mesh.voronoi(points, extent),
      vxs = [],
      vxIds = {},
      adjacencies = [],
      edges = [],
      triangles = [];
    for (let i = 0; i < voronoi.edges.length; i++) {
      let edge = voronoi.edges[i];
      if (edge == undefined) {
        continue;
      }
      let e0 = vxIds[edge[0]],
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
    let mesh = {
      adj: adjacencies,
      edges: edges,
      extent: extent,
      pts: points,
      tris: triangles,
      vor: voronoi,
      vxs: vxs
    };

    mesh.map = function (f) {
      let mapped = vxs.map(f);
      mapped.mesh = mesh;

      return mapped;
    };

    return mesh;
  }

  /**
   *
   * @param {landscapeObject} landscape
   * @param {function} f
   * @return {object} newH
   */
  static map(landscape, f) {
    let newH = landscape.map(f);
    newH.mesh = landscape.mesh;

    return newH;
  }

  /**
   *
   * @param {Array} segments
   * @return {Array} paths
   */
  static mergeSegments(segments) {
    let adj = {};
    for (let i = 0; i < segments.length; i++) {
      let seg = segments[i],
        a0 = adj[seg[0]] || [],
        a1 = adj[seg[1]] || [];
      a0.push(seg[1]);
      a1.push(seg[0]);
      adj[seg[0]] = a0;
      adj[seg[1]] = a1;
    }
    let done = [],
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
      let changed = false;
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
   * @param {meshObject} mesh
   * @param {number} number
   * @param {number} radius
   * @return {Array} newVals
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
        let m = mounts[j];
        newVals[i] += Math.pow(Math.exp(-((p[0] - m[0]) * (p[0] - m[0]) + (p[1] - m[1]) * (p[1] - m[1])) / (2 * radius * radius)), 2);
      }
    }

    return newVals;
  }

  /**
   *
   * @param {meshObject} mesh
   * @param {number} i
   * @return {Array} nbs
   */
  static neighbors(mesh, i) {
    let onbs = mesh.adj[i],
      nbs = [];
    for (let i = 0; i < onbs.length; i++) {
      nbs.push(onbs[i]);
    }

    return nbs;
  }

  /**
   * Rescale the heights to lie in the range 0-1
   * @param {landscapeObject} landscape
   * @return {*}
   */
  static normalize(landscape) {
    let lo = d3.min(landscape),
      hi = d3.max(landscape);

    return Mesh.map(landscape, function (x) {
      return (x - lo) / (hi - lo);
    });
  }

  /**
   * Normalize, then take the square root of the height value
   * to round off the tops of hills
   * @param {landscapeObject} landscape
   * @return {*}
   */
  static peaky(landscape) {
    return Mesh.map(Mesh.normalize(landscape), Math.sqrt);
  }

  /**
   *
   * @param {landscapeObject} landscape
   * @param {number} q
   * @return {number}
   */
  static quantile(landscape, q) {
    let sortedH = [];
    for (let i = 0; i < landscape.length; i++) {
      sortedH[i] = landscape[i];
    }
    sortedH.sort(d3.ascending);

    return d3.quantile(sortedH, q);
  }

  /**
   * Replace each height value with the average of its neighbours
   * to smooth the surface
   * @param {landscapeObject} landscape
   * @return {meshObject} newH
   */
  static relax(landscape) {
    let newH = Mesh.zero(landscape.mesh);
    for (let i = 0; i < landscape.length; i++) {
      let nbs = Mesh.neighbors(landscape.mesh, i);
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
   * @param {Array} path
   * @return {Array} newPath
   */
  static relaxPath(path) {
    let newPath = [path[0]];
    for (let i = 1; i < path.length - 1; i++) {
      let newPt = [
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
   * @param {landscapeObject} landscape
   * @param {number} q
   * @return {meshObject} newH
   */
  static setSeaLevel(landscape, q) {
    let newH = Mesh.zero(landscape.mesh),
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
   * @param {landscapeObject} landscape
   * @param {number} i
   * @return {*[]}
   */
  static triSlope(landscape, i) {
    let nbs = Mesh.neighbors(landscape.mesh, i);
    if (nbs.length != 3) {
      return [0, 0];
    }
    let p0 = landscape.mesh.vxs[nbs[0]],
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
   * @return {d3.voronoi} v
   */
  static voronoi(pts, extent) {
    let w = extent.width / 2,
      h = extent.height / 2;

    return d3.voronoi().extent([[-w, -h], [w, h]])(pts);
  }

  /**
   *
   * @param {meshObject} mesh
   * @return {landscapeObject} z
   */
  static zero(mesh) {
    let z = [];
    for (let i = 0; i < mesh.vxs.length; i++) {
      z[i] = 0;
    }
    z.mesh = mesh;

    return z;
  }

}

/**
 * @class Map
 */
class Map {

  /**
   *
   * @param {paramsObject} params
   * @param {ID3Selection} svg
   * @return {Map}
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
    this.viewErosionRate = false;

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
      newIslands: () => {
        this.clearMap();
        this.landscape = generateIslands(this.settings);
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

    return this;
  }

  /**
   *
   */
  addCity() {
    this.cities = this.cities || [];
    let score = this.cityScore(),
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
   * @return {object} score
   */
  cityScore() {
    let score = Mesh.map(Mesh.getFlux(this.landscape), Math.sqrt);
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
   * @return {Map}
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
    this.svg.select('g').selectAll('*').remove();

    return this;
  }

  /**
   * Generate all the coasts, borders, cities etc and draw them
   * @return {Map}
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
   * @param {string} cls
   * @param {Array} pts
   * @return {Map}
   */
  drawCircles(cls, pts) {
    let circles = this.svg.selectAll('circle.' + cls).data(pts);
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
   * @param {string} cls
   * @param {Array} pts
   * @return {Map}
   */
  drawLabels(cls, pts) {
    let cityTexts = this.svg.selectAll('text.city').data(this.cityLabels);
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


    let regionTexts = this.svg.selectAll('text.region').data(this.regionLabels);
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
   * @return {Map}
   */
  drawMap() {
    this.drawPaths('river', this.viewRivers ? this.rivers : []);
    this.drawPaths('coast', this.viewCoasts ? this.coasts : []);
    this.drawPaths('border', this.viewBorders ? this.borders : []);
    this.drawCircles('city', this.viewCities ? this.cities : []);
    this.drawLabels('city', this.viewLabels ? this.cityLabels : []);
    this.drawSlopes(this.viewSlopes ? this.landscape : Mesh.zero(this.landscape.mesh));
    // this.drawPoints(this.viewPoints ? this.landscape.mesh.pts : []);
    this.drawVoronoi(this.viewHeight ? this.landscape : Mesh.zero(this.landscape.mesh), -1, 1);

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
   * @param {string} cls
   * @param {Array} pts
   * @return {Map}
   */
  drawPaths(cls, pts) {
    let svgPaths = this.svg.selectAll('path.' + cls).data(pts);
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
   * @param {Array} pts
   * @return {Map}
   */
  drawPoints(pts) {
    let circle = this.svg.selectAll('circle').data(pts);
    circle.enter().append('circle');
    circle.exit().remove();
    this.svg.selectAll('circle')
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
    let score = this.cityScore();
    this.drawVoronoi(score, d3.max(score) - 0.5, d3.max(score) + 0.5);

    return this;
  }

  /**
   * Draw strokes which go up and right if the terrain slopes upwards from
   * left to right, and down and right if the terrain slopes downwards.
   * Similarly, the strokes on the 'near' side of hills should be longer than
   * those on the 'far' side.
   * @param {object} landscape
   * @return {Map}
   */
  drawSlopes(landscape) {
    let strokes = [],
      r = 0.25 / Math.sqrt(landscape.length),
      svgBase = this.svg.select('g');
    for (let i = 0; i < landscape.length; i++) {
      if (landscape[i] <= 0 || Mesh.isNearEdge(landscape.mesh, i)) {
        continue;
      }
      let nbs = Mesh.neighbors(landscape.mesh, i),
        s = 0,
        s2 = 0;
      nbs.push(i);
      for (let j = 0; j < nbs.length; j++) {
        let slopes = Mesh.triSlope(landscape, nbs[j]);
        s += slopes[0] / 10;
        s2 += slopes[1];
      }
      s /= nbs.length;
      s2 /= nbs.length;
      if (Math.abs(s) < runIf(0.1, 0.4)) {
        continue;
      }
      let l = r * runIf(1, 2) * (1 - 0.2 * Math.pow(Math.atan(s), 2)) * Math.exp(s2 / 100),
        x = landscape.mesh.vxs[i][0],
        y = landscape.mesh.vxs[i][1];
      if (Math.abs(l * s) > 2 * r) {
        let n = Math.floor(Math.abs(l * s / r));
        l /= n;
        if (n > 4) {
          n = 4;
        }
        for (let j = 0; j < n; j++) {
          let u = rNorm() * r,
            v = rNorm() * r;
          strokes.push([[x + u - l, y + v + l * s], [x + u + l, y + v - l * s]]);
        }
      } else {
        strokes.push([[x - l, y + l * s], [x + l, y - l * s]]);
      }
    }

    let lines = this.svg.selectAll('line.slope').data(strokes);
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
   * @param {object} landscape
   * @param {number|null} lo
   * @param {number|null} hi
   */
  drawVoronoi(landscape, lo = null, hi = null) {
    if (hi === null) {
      hi = d3.max(landscape) + 1e-9;
    }
    if (lo === null) {
      lo = d3.min(landscape) - 1e-9;
    }
    let mappedVals = landscape.map(function (x) {
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
   * @param {object} landscape
   * @param {object} territories
   * @return {Map}
   */
  getBorders(landscape, territories) {
    let edges = [];
    for (let i = 0; i < territories.mesh.edges.length; i++) {
      let e = territories.mesh.edges[i];
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
   * @return {Map}
   */
  getCities() {
    let n = this.settings.numCities;
    for (let i = 0; i < n; i++) {
      let score = this.cityScore(),
        newCity = d3.scan(score, d3.descending);
      this.cities.push(newCity);
    }

    return this;
  }

  /**
   * Generate and set the coasts for the landscape
   * @param {object} landscape
   * @return {Map}
   */
  getCoasts(landscape) {
    this.coasts = Mesh.contour(landscape, 0);

    return this;
  }

  /**
   * Generate and set the labels for the landscape
   * @param {object} landscape
   * @return {Map}
   */
  getLabels(landscape) {
    let numTer = this.settings.numTerritories,
      avoids = [this.rivers, this.coasts, this.borders],
      lang = makeRandomLanguage();
    this.cityLabels = [];

    let penalty = (label) => {
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
        town = i >= numTer,
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
   * @param {object} landscape
   * @param {number} limit
   * @return {Map}
   */
  getRivers(landscape, limit = 0.01) {
    let dh = Mesh.downhill(landscape),
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
   * @return {Map}
   */
  getTerritories(landscape) {
    this.territories = [];
    let n = this.settings.numTerritories;
    if (n > this.cities.length) {
      n = this.cities.length;
    }
    let flux = Mesh.getFlux(landscape),
      queue = new PriorityQueue({
        comparator: function (a, b) {
          return a.score - b.score;
        }
      });

    let weight = (u, v) => {
      let horiz = Mesh.distance(landscape.mesh, u, v),
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
      let nbs = Mesh.neighbors(landscape.mesh, u.vx);
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
   * @param {object} terr
   * @param {object} city
   * @param {boolean} landOnly
   * @return {*[]}
   */
  getTerritoryCenter(terr, city, landOnly) {
    let x = 0, y = 0, n = 0;
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
   * @param {Array} path
   * @return {string}
   */
  static makeD3Path(path) {
    let p = d3.path();
    p.moveTo(1000 * path[0][0], 1000 * path[0][1]);
    for (let i = 1; i < path.length; i++) {
      p.lineTo(1000 * path[i][0], 1000 * path[i][1]);
    }

    return p.toString();
  }

}
