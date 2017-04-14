
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
    let x = col,
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
    let x = Math.trunc(Math.round(hex.x)),
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
}
