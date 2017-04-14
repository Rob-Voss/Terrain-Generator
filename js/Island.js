// Park-Miller-Carta Pseudo-Random Number Generator
function PRNG() {
  this.seed = 1;
  this.next = () => this.gen() / 2147483647;
  this.nextRange = (min, max) => min + ((max - min) * this.next());
  this.gen = () => this.seed = (this.seed * 16807) % 2147483647;
}

/**
 *
 * @param {HTMLCanvasElement} canvas
 * @param {number} baseX
 * @param {number} baseY
 * @param {number} seed
 * @returns {ImageData}
 */
function perlinNoise(canvas, baseX, baseY, seed) {
  let rand = new PRNG(),
    ctx = canvas.getContext('2d'),
    imageData = ctx.createImageData(canvas.width, canvas.height),
    data = imageData.data,
    simplexR = new SimplexNoise(rand),
    simplexG = new SimplexNoise(rand),
    simplexB = new SimplexNoise(rand);

  simplexR.setSeed(seed);
  simplexG.setSeed(seed + 1);
  simplexB.setSeed(seed + 2);

  let pos, cr, cg, cb, gray;
  for (let y = 0; y < canvas.height; y++) {
    for (let x = 0; x < canvas.width; x++) {
      pos = (x + y * canvas.width) * 4;

      cr = Math.floor(((simplexR.noise(x / baseX, y / baseY) + 1) * 0.5) * 255);
      cg = Math.floor(((simplexG.noise(x / baseX, y / baseY) + 1) * 0.5) * 255);
      cb = Math.floor(((simplexB.noise(x / baseX, y / baseY) + 1) * 0.5) * 255);

      gray = (cr + cg + cb) / 3;

      data[pos + 0] = gray;
      data[pos + 1] = gray;
      data[pos + 2] = gray;
      data[pos + 3] = 255;
    }
  }
  ctx.putImageData(imageData, 0, 0);

  return imageData;
}

/**
 * @class Island
 *
 * @property {Paper.Layer} cellsLayer
 * @property {boolean} debug
 * @property {Paper.Layer} debugLayer
 * @property {Voronoi.Diagram} diagram
 * @property {ImageData} perlin
 * @property {Paper.Layer} riversLayer
 * @property {number} seed
 * @property {Array} sites
 * @property {Voronoi} voronoi
 */
class Island {

  /**
   *
   *
   */
  constructor() {
    this.debug = false;
    this.diagram = null;
    this.perlin = null;
    this.sites = [];
    this.voronoi = new Voronoi();
  }

  /**
   *
   */
  assignBiomes() {
    for (let i = 0; i < this.diagram.cells.length; i++) {
      this.diagram.cells[i].biome = Island.getBiome(this.diagram.cells[i]);
    }
  }

  /**
   * Calculate moisture.
   * Freshwater sources spread moisture: rivers and lakes (not ocean).
   */
  assignMoisture() {
    let queue = [];
    // Lake and river
    for (let i = 0; i < this.diagram.cells.length; i++) {
      let cell = this.diagram.cells[i];
      if ((cell.water || cell.river) && !cell.ocean) {
        cell.moisture = (cell.water ? 1 : 0.9);
        if (!cell.ocean) {
          queue.push(cell);
        }
      }
    }

    while (queue.length > 0) {
      let cell = queue.shift(),
        neighbors = cell.getNeighborIds();
      for (let i = 0; i < neighbors.length; i++) {
        let nId = neighbors[i],
          neighbor = this.diagram.cells[nId],
          newMoisture = cell.moisture * 0.9;
        if (neighbor.moisture == null || newMoisture > neighbor.moisture) {
          neighbor.moisture = newMoisture;
          queue.push(neighbor);
        }
      }
    }

    // Ocean
    for (let i = 0; i < this.diagram.cells.length; i++) {
      let cell = this.diagram.cells[i];
      if (cell.ocean) {
        cell.moisture = 1;
      }
    }
  }

  /**
   *
   */
  assignOceanCoastAndLand() {
    // Water
    let queue = [];
    for (let i = 0; i < this.diagram.cells.length; i++) {
      let cell = this.diagram.cells[i];
      cell.elevation = this.getElevation(cell.site);
      cell.water = (cell.elevation <= 0);
      let numWater = 0;
      for (let j = 0; j < cell.halfedges.length; j++) {
        let hedge = cell.halfedges[j];
        // Border
        if (hedge.edge.rSite == null) {
          cell.border = true;
          cell.ocean = true;
          cell.water = true;
          if (cell.elevation > numWater) {
            cell.elevation = numWater;
          }
          queue.push(cell);
        }
      }
    }

    // Ocean
    while (queue.length > 0) {
      let cell = queue.shift(),
        neighbors = cell.getNeighborIds();
      for (let i = 0; i < neighbors.length; i++) {
        let nId = neighbors[i],
          neighbor = this.diagram.cells[nId];
        if (neighbor.water && !neighbor.ocean) {
          neighbor.ocean = true;
          queue.push(neighbor);
        }
      }
    }

    // Coast
    for (let i = 0; i < this.diagram.cells.length; i++) {
      let cell = this.diagram.cells[i],
        numOcean = 0,
        neighbors = cell.getNeighborIds();
      for (let j = 0; j < neighbors.length; j++) {
        let nId = neighbors[j],
          neighbor = this.diagram.cells[nId];
        if (neighbor.ocean) {
          numOcean++;
        }
      }
      cell.coast = (numOcean > 0) && (!cell.water);
      cell.beach = (cell.coast && cell.elevation < this.cliffsThreshold);
    }

    // Cliff
    for (let i = 0; i < this.diagram.edges.length; i++) {
      let edge = this.diagram.edges[i];
      if (edge.lSite != null && edge.rSite != null) {
        let lCell = this.diagram.cells[edge.lSite.voronoiId],
          rCell = this.diagram.cells[edge.rSite.voronoiId],
          elevationDiff = Math.abs(Island.getRealElevation(lCell) - Island.getRealElevation(rCell));

        edge.cliff = (!(lCell.water && rCell.water) && (elevationDiff >= this.cliffsThreshold));
      }
    }
  }

  /**
   *
   */
  assignRivers() {
    for (let i = 0; i < this.numRivers;) {
      let cell = this.diagram.cells[Island.getRandomInt(0, this.diagram.cells.length - 1)];
      if (!cell.coast) {
        if (this.setAsRiver(cell, 1)) {
          cell.source = true;
          i++;
        }
      }
    }
  }

  /**
   *
   * @param {Voronoi.cell} cell
   * @returns {number}
   */
  static cellArea(cell) {
    let area = 0,
      halfEdges = cell.halfedges,
      iHalfEdge = halfEdges.length,
      halfEdge,
      p1, p2;
    while (iHalfEdge--) {
      halfEdge = halfEdges[iHalfEdge];
      p1 = halfEdge.getStartpoint();
      p2 = halfEdge.getEndpoint();
      area += p1.x * p2.y;
      area -= p1.y * p2.x;
    }
    area /= 2;

    return area;
  }

  /**
   *
   * @param {Voronoi.cell} cell
   * @returns {{x: number, y: number}}
   */
  static cellCentroid(cell) {
    let x = 0,
      y = 0,
      halfEdges = cell.halfedges,
      iHalfEdge = halfEdges.length,
      halfEdge,
      v, p1, p2, pt;
    while (iHalfEdge--) {
      halfEdge = halfEdges[iHalfEdge];
      p1 = halfEdge.getStartpoint();
      p2 = halfEdge.getEndpoint();
      v = p1.x * p2.y - p2.x * p1.y;
      x += (p1.x + p2.x) * v;
      y += (p1.y + p2.y) * v;
    }
    v = Island.cellArea(cell) * 6;
    pt = new Point(x / v, y / v);

    return pt;
  }

  /**
   *
   * @param {Array} sites
   */
  compute(sites) {
    this.sites = sites;
    this.voronoi.recycle(this.diagram);
    let bBox = {
      xl: 0,
      xr: this.width,
      yt: 0,
      yb: this.height
    };
    this.diagram = this.voronoi.compute(sites, bBox);
  }

  /**
   *
   * @param {Point} a
   * @param {Point} b
   * @returns {number}
   */
  static distance(a, b) {
    let dx = a.x - b.x,
      dy = a.y - b.y;

    return Math.sqrt(dx * dx + dy * dy);
  }

  /**
   *
   * @param {Voronoi.cell} cell
   */
  fillLake(cell) {
    // If the lake has an exit river it can not longer be filled
    if (cell.exitRiver == null) {
      let exitRiver = null,
        exitSource = null,
        lake = [],
        queue = [];
      queue.push(cell);

      while (queue.length > 0) {
        let c = queue.shift();
        lake.push(c);
        let neighbors = c.getNeighborIds();
        for (let i = 0; i < neighbors.length; i++) {
          let nId = neighbors[i],
            neighbor = this.diagram.cells[nId];
          if (neighbor.water && !neighbor.ocean) {
            // Water cell from the same lake
            if (neighbor.lakeElevation == null || neighbor.lakeElevation < c.lakeElevation) {
              neighbor.lakeElevation = c.lakeElevation;
              queue.push(neighbor);
            }
          } else {
            // Ground cell adjacent to the lake
            if (c.elevation < neighbor.elevation) {
              if (neighbor.elevation - c.lakeElevation < 0) {
                // We fill the ground with water
                neighbor.water = true;
                neighbor.lakeElevation = c.lakeElevation;
                queue.push(neighbor);
              }
            } else {
              // neighbor.source = true;
              // We found a new exit for the lake
              if (exitRiver == null || exitRiver.elevation > neighbor.elevation) {
                exitSource = c;
                exitRiver = neighbor;
              }
            }
          }
        }
      }

      if (exitRiver != null) {
        // We start the exit river
        exitSource.river = true;
        exitSource.nextRiver = exitRiver;
        this.setAsRiver(exitRiver, 2);
        // We mark all the lake as having an exit river
        while (lake.length > 0) {
          let c = lake.shift();
          c.exitRiver = exitRiver;
        }
      }
    }
  }

  /**
   *
   * @param {Voronoi.cell} cell
   * @returns {string}
   */
  static getBiome(cell) {
    if (cell.ocean) {
      return 'OCEAN';
    } else if (cell.water) {
      if (Island.getRealElevation(cell) < 0.05) {
        return 'MARSH';
      }
      if (Island.getRealElevation(cell) > 0.4) {
        return 'ICE';
      }
      return 'LAKE';
    } else if (cell.beach) {
      return 'BEACH';
    } else if (cell.elevation > 0.4) {
      if (cell.moisture > 0.50) {
        return 'SNOW';
      } else if (cell.moisture > 0.33) {
        return 'TUNDRA';
      } else if (cell.moisture > 0.16) {
        return 'BARE';
      } else {
        return 'SCORCHED';
      }
    } else if (cell.elevation > 0.3) {
      if (cell.moisture > 0.66) {
        return 'TAIGA';
      } else if (cell.moisture > 0.33) {
        return 'SHRUBLAND';
      } else {
        return 'TEMPERATE_DESERT';
      }
    } else if (cell.elevation > 0.15) {
      if (cell.moisture > 0.83) {
        return 'TEMPERATE_RAIN_FOREST';
      } else if (cell.moisture > 0.50) {
        return 'TEMPERATE_DECIDUOUS_FOREST';
      } else if (cell.moisture > 0.16) {
        return 'GRASSLAND';
      } else {
        return 'TEMPERATE_DESERT';
      }
    } else {
      if (cell.moisture > 0.66) {
        return 'TROPICAL_RAIN_FOREST';
      } else if (cell.moisture > 0.33) {
        return 'TROPICAL_SEASONAL_FOREST';
      } else if (cell.moisture > 0.16) {
        return 'GRASSLAND';
      } else {
        return 'SUBTROPICAL_DESERT';
      }
    }
  }

  /**
   *
   * @param {Voronoi.cell} cell
   */
  getCellColor(cell) {
    let c = Island.DISPLAY_COLORS[cell.biome].clone();
    c.brightness = c.brightness - this.getShade(cell);

    return c;
  }

  /**
   * The Perlin-based island combines perlin noise with the radius
   *
   * @param {Point} point
   * @returns {number}
   */
  getElevation(point) {
    let x = 2 * (point.x / this.width - 0.5),
      y = 2 * (point.y / this.height - 0.5),
      distance = Math.sqrt(x * x + y * y),
      c = this.getPerlinValue(point),
      // Multiple small islands and one big one
      elevation1 = c - (0.3 + 0.3 * distance * distance),
      // One island
      elevation = c - distance;

    return elevation1;
    // return elevation;
  }

  /**
   *
   * @param {Point} point
   * @returns {number}
   */
  getPerlinValue(point) {
    let x = ((point.x / this.width) * this.perlin.width) | 0,
      y = ((point.y / this.height) * this.perlin.height) | 0,
      pos = (x + y * this.perlin.width) * 4,
      data = this.perlin.data,
      val = data[pos] << 16 | data[pos + 1] << 8 | data[pos + 2]; // rgb to hex

    return (val & 0xff) / 255.0;
  }

  /**
   *
   * @param {number} min
   * @param {number} max
   * @returns {number}
   */
  static getRandomInt(min, max) {
    return Math.floor(Math.random() * (max - min + 1)) + min;
  }

  /**
   *
   * @param {Voronoi.cell} cell
   * @returns {number}
   */
  static getRealElevation(cell) {
    if (cell.water && cell.lakeElevation != null) {
      return cell.lakeElevation;
    } else if (cell.water && cell.elevation < 0) {
      return 0;
    } else {
      return cell.elevation;
    }
  }

  /**
   *
   * @param {Voronoi.cell} cell
   * @returns {number}
   */
  getShade(cell) {
    if (this.shading == 0) {
      return 0;
    } else if (cell.ocean) {
      return (this.shadeOcean ? -cell.elevation : 0);
    } else if (cell.water) {
      return 0;
    } else {
      let lowerCell = null,
        upperCell = null,
        neighbors = cell.getNeighborIds();
      for (let j = 0; j < neighbors.length; j++) {
        let nId = neighbors[j],
          neighbor = this.diagram.cells[nId];
        if (lowerCell == null || neighbor.elevation < lowerCell.elevation) {
          lowerCell = neighbor;
        }
        if (upperCell == null || neighbor.elevation > upperCell.elevation) {
          upperCell = neighbor;
        }
      }

      let angleRadian = Math.atan2(upperCell.site.x - lowerCell.site.x, upperCell.site.y - lowerCell.site.y),
        angleDegree = angleRadian * (180 / Math.PI),
        diffElevation = (Island.getRealElevation(upperCell) - Island.getRealElevation(lowerCell));

      if (diffElevation + this.shading < 1) {
        diffElevation = diffElevation + this.shading;
      }

      return ((Math.abs(angleDegree) / 180) * diffElevation);
    }
  }

  /**
   *
   * @param {HTMLCanvasElement} islandCanvas
   * @param {HTMLCanvasElement} perlinCanvas
   * @param {Island.defaultUserConfig} userConfig
   * @returns {Island}
   */
  init(islandCanvas, perlinCanvas, userConfig = {}) {
    this.width = userConfig.width || Island.defaultUserConfig.width;
    this.height = userConfig.height || Island.defaultUserConfig.height;
    this.perlinWidth = userConfig.perlinWidth || (this.width / 3);
    this.perlinHeight = userConfig.perlinHeight || (this.height / 3);
    this.allowDebug = userConfig.allowDebug || false;
    this.numSites = userConfig.numSites || ((this.width * this.height) / 100);
    this.sitesDistribution = userConfig.sitesDistribution || 'hexagon';
    this.sitesRandomisation = userConfig.sitesRandomisation || 80;
    this.numGraphRelaxation = userConfig.numGraphRelaxation || 0;
    this.cliffsThreshold = userConfig.cliffsThreshold || 0.15;
    this.lakesThreshold = userConfig.lakesThreshold || 0.005;
    this.numRivers = userConfig.numRivers || (this.numSites / 200);
    this.maxRiversSize = userConfig.maxRiversSize || 4;
    this.shading = userConfig.shading || 0.35;
    this.shadeOcean = userConfig.shadeOcean || true;
    this.cellsLayer = new paper.Layer({name: 'cell'});
    this.riversLayer = new paper.Layer({name: 'rivers'});
    this.debugLayer = new paper.Layer({name: 'debug', visible: false});
    this.seed = Math.random();

    this.islandCanvas = islandCanvas;
    this.islandCanvas.width = this.width;
    this.islandCanvas.height = this.height;
    this.perlinCanvas = perlinCanvas;
    this.perlinCanvas.width = this.perlinWidth;
    this.perlinCanvas.height = this.perlinHeight;
    this.perlin = perlinNoise(this.perlinCanvas, 64, 64, this.seed);
    this.randomSites();

    this.assignOceanCoastAndLand();
    this.assignRivers();
    this.assignMoisture();
    this.assignBiomes();

    this.render();

    return this;
  }

  /**
   *
   */
  randomSites() {
    let sites = [];

    // Create vertices
    if (this.sitesDistribution == 'random') {
      for (let i = 0; i < this.numSites; i++) {
        sites.push(new Point(
          Math.round(Math.random() * this.width),
          Math.round(Math.random() * this.height)
        ));
      }
    } else {
      let delta = Math.sqrt(this.width * this.height / this.numSites),
        rand = this.sitesRandomisation * delta / 100,
        x = 0,
        y = 0;
      for (let i = 0; i < this.numSites; i++) {
        sites.push(new Point(
          Math.max(Math.min(Math.round(x * delta + (Math.random() * rand)), this.width), 0),
          Math.max(Math.min(Math.round(y * delta + (Math.random() * rand)), this.height), 0)
        ));
        x = x + 1;
        if (x * delta > this.width) {
          x = (y % 2 == 1 || this.sitesDistribution == 'square' ? 0 : 0.5);
          y = y + 1;
        }
      }
    }

    this.compute(sites);

    for (let i = 0; i < this.numGraphRelaxation; i++) {
      this.relaxSites();
    }
  }

  /**
   *
   */
  relaxSites() {
    if (!this.diagram) {
      return;
    }
    let cells = this.diagram.cells,
      iCell = cells.length,
      cell,
      site, sites = [],
      rn, dist,
      p = 1 / iCell * 0.1;
    while (iCell--) {
      cell = cells[iCell];
      rn = Math.random();
      // Probability of apoptosis
      // (the death of cells that occurs as a normal and controlled part of an organism's growth or development.)
      if (rn < p) {
        continue;
      }
      site = Island.cellCentroid(cell);
      dist = Island.distance(site, cell.site);
      // Don't relax too fast
      if (dist > 2) {
        site.x = (site.x + cell.site.x) / 2;
        site.y = (site.y + cell.site.y) / 2;
      }
      // Probability of mitosis
      if (rn > (1 - p)) {
        dist /= 2;
        sites.push(new Point(
          site.x + (site.x - cell.site.x) / dist,
          site.y + (site.y - cell.site.y) / dist
        ));
      }
      sites.push(site);
    }
    this.compute(sites);
  }

  /**
   *
   */
  render() {
    if (!this.diagram) {
      return;
    }

    this.renderCells();
    this.renderRivers();
    this.renderEdges();
    this.renderSites();

    paper.view.draw();
  }

  /**
   *
   */
  renderCells() {
    this.cellsLayer.activate();
    for (let cellId in this.diagram.cells) {
      let cell = this.diagram.cells[cellId],
        color = this.getCellColor(cell),
        cellPath = new Path();
      cellPath.strokeWidth = 1;
      cellPath.strokeColor = color;
      cellPath.fillColor = color;
      let start = cell.halfedges[0].getStartpoint(),
        startPt = new Point(start.x, start.y);
      cellPath.add(startPt);
      for (let i = 0; i < cell.halfedges.length; i++) {
        let halfEdge = cell.halfedges[i],
          end = halfEdge.getEndpoint(),
          endPt = new Point(end.x, end.y);
        cellPath.add(endPt);
      }
      cellPath.closed = true;
    }
  }

  /**
   *
   */
  renderEdges() {
    if (this.allowDebug) {
      this.debugLayer.activate();
      let edges = this.diagram.edges,
        iEdge = edges.length,
        edge, v;
      while (iEdge--) {
        edge = edges[iEdge];
        let edgePath = new Path();
        edgePath.strokeWidth = 1;

        if (edge.cliff) {
          edgePath.strokeWidth = 1;
          edgePath.strokeCap = 'round';
          edgePath.strokeColor = Island.DISPLAY_COLORS.ROCK;
        } else {
          edgePath.strokeWidth = 1;
          edgePath.strokeColor = '#000';
        }
        v = edge.va;
        edgePath.add(new Point(v.x, v.y));
        v = edge.vb;
        edgePath.add(new Point(v.x, v.y));
      }
    }
  }

  /**
   *
   */
  renderRivers() {
    for (let cellId in this.diagram.cells) {
      let cell = this.diagram.cells[cellId];
      if (cell.nextRiver) {
        this.riversLayer.activate();
        let riverPath = new Path();
        riverPath.strokeWidth = Math.min(cell.riverSize, this.maxRiversSize);
        let riverColor = Island.DISPLAY_COLORS.RIVER.clone();
        riverColor.brightness = riverColor.brightness - this.getShade(cell);
        riverPath.strokeColor = riverColor;
        riverPath.strokeCap = 'round';
        if (cell.water) {
          riverPath.add(new Point(
            cell.site.x + (cell.nextRiver.site.x - cell.site.x) / 2,
            cell.site.y + (cell.nextRiver.site.y - cell.site.y) / 2
          ));
        } else {
          riverPath.add(new Point(
            cell.site.x,
            cell.site.y
          ));
        }
        if (cell.nextRiver && !cell.nextRiver.water) {
          riverPath.add(new Point(
            cell.nextRiver.site.x,
            cell.nextRiver.site.y
          ));
        } else {
          riverPath.add(new Point(
            cell.site.x + (cell.nextRiver.site.x - cell.site.x) / 2,
            cell.site.y + (cell.nextRiver.site.y - cell.site.y) / 2
          ));
        }
      }
      // source :
      if (this.allowDebug && cell.source) {
        this.debugLayer.activate();
        let circle = new Path.Circle(new Point(
          cell.site.x,
          cell.site.y
        ), 3);
        circle.fillColor = Island.DISPLAY_COLORS.SOURCE;
      }
    }
  }

  /**
   *
   */
  renderSites() {
    if (this.allowDebug) {
      this.debugLayer.activate();
      // Sites
      let sites = this.sites,
        iSite = sites.length;
      while (iSite--) {
        let v = sites[iSite],
          circle = new Path.Circle(new Point(v.x, v.y), 1);
        circle.fillColor = '#0f0';
      }

      // Values
      for (let i = 0; i < this.diagram.cells.length; i++) {
        let cell = this.diagram.cells[i],
          text = new PointText(new Point(cell.site.x, cell.site.y));
        text.fillColor = '#f00';
        text.fontSize = '8px';
        text.content = Math.ceil(Island.getRealElevation(cell) * 100) + "\r\n" + Island.getBiome(cell);
      }
    }
  }

  /**
   *
   * @param {Voronoi.cell} cell
   * @param {number} size
   * @returns {boolean}
   */
  setAsRiver(cell, size) {
    if (!cell.water && !cell.river) {
      cell.river = true;
      cell.riverSize = size;
      let lowerCell = null,
        neighbors = cell.getNeighborIds();
      // We choose the lowest neighbour cell
      for (let j = 0; j < neighbors.length; j++) {
        let nId = neighbors[j],
          neighbor = this.diagram.cells[nId];
        if (lowerCell == null || neighbor.elevation < lowerCell.elevation) {
          lowerCell = neighbor;
        }
      }
      if (lowerCell.elevation < cell.elevation) {
        // We continue the river to the next lowest cell
        this.setAsRiver(lowerCell, size);
        cell.nextRiver = lowerCell;
      } else {
        // We are in a hole, so we create a lake
        cell.water = true;
        this.fillLake(cell);
      }
    } else if (cell.water && !cell.ocean) {
      // We ended in a lake, the water level rise
      cell.lakeElevation = Island.getRealElevation(cell) + (this.lakesThreshold * size);
      this.fillLake(cell);
    } else if (cell.river) {
      // We ended in another river, the river size increase
      cell.riverSize++;
      let nextRiver = cell.nextRiver;
      while (nextRiver) {
        nextRiver.riverSize++;
        nextRiver = nextRiver.nextRiver;
      }
    }

    return cell.river;
  }

  /**
   *
   */
  toggleDebug() {
    this.debug = !this.debug;
    this.debugLayer.visible = this.debug;
  }
}

/**
 *
 * @typedef {Object} defaultUserConfig
 * @property {number} width
 * @property {number} height
 * @property {number} perlinWidth
 * @property {number} perlinHeight
 * @property {boolean} allowDebug
 * @property {number} numSites
 * @property {string} sitesDistribution
 * @property {number} sitesRandomisation
 * @property {number} numGraphRelaxation
 * @property {number} cliffsThreshold
 * @property {number} lakesThreshold
 * @property {number} numRivers
 * @property {number} maxRiversSize
 * @property {number} shading
 * @property {boolean} shadeOcean
 */
Island.defaultUserConfig = {
  width: 500,
  height: 500,
  perlinWidth: 256,
  perlinHeight: 256,
  // If set to true, you can click on the map to enter "debug" mode.
  // Warning : debug mode is slow to initialize, set to false for faster rendering.
  allowDebug: false,
  // Number of voronoi cells
  numSites: 10000,
  // Distribution of the site : random, square or hexagon
  sitesDistribution: 'hexagon',
  // Will move each site in a random way (in %), for the square or hexagon distribution to look more random
  sitesRandomisation: 80,
  // Number of time we apply the relaxation algo to the voronoi graph (slow !), for the random distribution to look less random
  numGraphRelaxation: 0,
  cliffsThreshold: 0.15,
  // lake elevation will increase by this value (* the river size) when a new river end inside
  lakesThreshold: 0.005,
  numRivers: (10000 / 200),
  maxRiversSize: 4,
  shading: 0.35,
  shadeOcean: true
};

/**
 *
 * @typedef {Object} DISPLAY_COLORS
 * @property {Color} OCEAN
 * @property {Color} BEACH
 * @property {Color} LAKE
 * @property {Color} RIVER
 * @property {Color} SOURCE
 * @property {Color} MARSH
 * @property {Color} ICE
 * @property {Color} ROCK
 * @property {Color} LAVA
 * @property {Color} SNOW
 * @property {Color} TUNDRA
 * @property {Color} BARE
 * @property {Color} SCORCHED
 * @property {Color} TAIGA
 * @property {Color} SHRUBLAND
 * @property {Color} TEMPERATE_DESERT
 * @property {Color} TEMPERATE_RAIN_FOREST
 * @property {Color} TEMPERATE_DECIDUOUS_FOREST
 * @property {Color} GRASSLAND
 * @property {Color} TROPICAL_RAIN_FOREST
 * @property {Color} TROPICAL_SEASONAL_FOREST
 * @property {Color} SUBTROPICAL_DESERT
 */
Island.DISPLAY_COLORS = {
  OCEAN: new paper.Color('#82caff'),
  BEACH: new paper.Color('#ffe98d'),
  LAKE: new paper.Color('#2f9ceb'),
  RIVER: new paper.Color('#369eea'),
  SOURCE: new paper.Color('#00f'),
  MARSH: new paper.Color('#2ac6d3'),
  ICE: new paper.Color('#b3deff'),
  ROCK: new paper.Color('#535353'),
  LAVA: new paper.Color('#e22222'),
  SNOW: new paper.Color('#f8f8f8'),
  TUNDRA: new paper.Color('#ddddbb'),
  BARE: new paper.Color('#bbbbbb'),
  SCORCHED: new paper.Color('#999999'),
  TAIGA: new paper.Color('#ccd4bb'),
  SHRUBLAND: new paper.Color('#c4ccbb'),
  TEMPERATE_DESERT: new paper.Color('#e4e8ca'),
  TEMPERATE_RAIN_FOREST: new paper.Color('#a4c4a8'),
  TEMPERATE_DECIDUOUS_FOREST: new paper.Color('#b4c9a9'),
  GRASSLAND: new paper.Color('#c4d4aa'),
  TROPICAL_RAIN_FOREST: new paper.Color('#9cbba9'),
  TROPICAL_SEASONAL_FOREST: new paper.Color('#a9cca4'),
  SUBTROPICAL_DESERT: new paper.Color('#e9ddc7')
};
