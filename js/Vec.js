const toRadian = Math.PI / 180,
  toDegree = 180 / Math.PI;

class Point {

  /**
   * Simple Point
   * @constructor
   *
   * @param {number} x
   * @param {number} y
   * @return {Point}
   */
  constructor(x = 0, y = 0) {
    this.x = x;
    this.y = y;

    return this;
  }
}

class Vec {

  /**
   * A 2D vector utility
   * @constructor
   *
   * @param {number} x
   * @param {number} y
   * @param {number} vx
   * @param {number} vy
   * @param {number} ax
   * @param {number} ay
   * @return {Vec}
   */
  constructor(x = 0, y = 0, vx = 0, vy = 0, ax = 0, ay = 0) {
    this.x = x;
    this.y = y;
    this.vx = vx;
    this.vy = vy;
    this.ax = ax;
    this.ay = ay;
    this.angle = Math.atan2(y, x);
    this.direction = 0;

    return this;
  }

  /**
   * Adds a vector to this one
   * @param {Vec} v The vector to add to this one
   * @return {Vec} Returns itself.
   */
  addVecTo(v) {
    return new Vec(this.x + v.x, this.y + v.y, this.vx + v.vx, this.vy + v.vy, this.ax + v.ax, this.ay + v.ay);
  }

  /**
   * Adds two vectors to each other and stores the result in this vector
   * @param {Vec} a
   * @param {Vec} b
   * @return {Vec} Returns itself.
   */
  addVecsTo(a, b) {
    this.x = a.x + b.x;
    this.y = a.y + b.y;

    return this;
  }

  /**
   * Adds a scalar value to the x and y components of this vector
   * @param {number} s The scalar value to add
   * @return {Vec} Returns itself.
   */
  addScalar(s) {
    this.x += s;
    this.y += s;

    return this;
  }

  /**
   * This will add the velocity x,y to the position x,y
   * @return {Vec}
   */
  advance(speed) {
    let oldPos = this.clone();
    this.x += this.vx;
    this.y += this.vy;
    this.angle = Math.atan2(this.y, this.x);
    this.ax = speed * Math.cos(this.angle);
    this.ay = speed * Math.sin(this.angle);
    this.direction = oldPos.angleBetween(this);

    return this;
  }

  /**
   * Get the angle of this Vec
   *       90
   *       ^
   * 180 <-|-> 0/360
   *       v
   *       270
   * @param {boolean} inDegree
   * @return {number}
   */
  getAngle(inDegree = false) {
    this.angle = Math.atan2(this.y, this.x);

    return (inDegree) ? this.angle * toDegree : this.angle;
  }

  /**
   * Calculate angle between any two vectors.
   * @param {Vec} v First vec
   * @param {boolean} inDegree
   * @return {number} Angle between vectors.
   */
  angleBetween(v, inDegree = false) {
    let v1 = this.clone(),
      v2 = v.clone(),
      angle = Math.atan2(v2.subByVec(v1).y, v2.subByVec(v1).x);

    return (inDegree) ? angle * toDegree : angle;
  }

  /**
   * Ceils the vector components
   * @return {Vec} Returns itself.
   */
  ceil() {
    this.x = Math.ceil(this.x);
    this.y = Math.ceil(this.y);

    return this;
  }

  /**
   * Clamps the vectors components to be between min and max
   * This function assumes min < max, if this assumption
   * isn't true it will not operate correctly
   * @param {Vec} min The minimum value a component can be
   * @param {Vec} max The maximum value a component can be
   * @return {Vec} Returns itself.
   */
  clamp(min, max) {
    if (this.x < min.x) {
      this.x = min.x;
    } else if (this.x > max.x) {
      this.x = max.x;
    }

    if (this.y < min.y) {
      this.y = min.y;
    } else if (this.y > max.y) {
      this.y = max.y;
    }

    return this;
  }

  /**
   * Creates a new instance of Vector, with the same components as this vector
   * @return {Vec} Returns a new Vector with the same values
   */
  clone() {
    return new Vec(this.x, this.y, this.vx, this.vy, this.ax, this.ay);
  }

  /**
   * Copies the passed vector's components to this vector
   * @param {Vec} v The vector to copy the values from
   * @return {Vec} Returns itself.
   */
  copy(v) {
    this.x = v.x;
    this.y = v.y;

    return this;
  }

  /**
   * Calculate the cross product of this and another vector.
   * @param {Vec} v A vector
   * @return {number} The cross product
   */
  crossProd(v) {
    return this.x * v.y - this.y * v.x;
  }

  /**
   * Calculates the square distance to the passed vector
   * @param {Vec} v The vector to check distance to
   * @return {number} The square distance
   */
  distanceToSquared(v) {
    let dx = this.x - v.x,
      dy = this.y - v.y;

    return dx * dx + dy * dy;
  }

  /**
   * Calculates the distance to the passed vector
   * @param {Vec} v The vector to check distance to
   * @return {number} The distance
   */
  distanceTo(v) {
    return Math.sqrt(this.distanceToSquared(v));
  }

  /**
   * Divides the x and y components of this vector by a scalar value
   * @param {number} s The value to divide by
   * @return {Vec} Returns itself.
   */
  divideBy(s) {
    if (s !== 0) {
      this.x /= s;
      this.y /= s;
    } else {
      this.x = 0;
      this.y = 0;
    }

    return this;
  }

  /**
   * Performs the dot product between this vector and
   * the passed one and returns the result
   * @param {Vec} v
   * @return {number} Returns the dot product
   */
  dot(v) {
    return (this.x * v.x) + (this.y * v.y);
  }

  /**
   * Checks if this vector is equal to another
   * @param {Vec} v The vector to compare with
   * @return {boolean}
   */
  equals(v) {
    return ((v.x === this.x) && (v.y === this.y));
  }

  /**
   * Floors the vector components
   * @return {Vec} Returns itself.
   */
  floor() {
    this.x = Math.floor(this.x);
    this.y = Math.floor(this.y);

    return this;
  }

  /**
   * Get a point at a % point between this Vec and another
   * @param {Vec} v
   * @param {number} p
   * @return {Vec} .
   */
  getPointBetween(v, p) {
    let blend = p / 100,
      x = this.x + blend * (v.x - this.x),
      y = this.y + blend * (v.y - this.y);

    return new Vec(x, y);
  }

  /**
   * Calculates the square length of the vector
   * @return {number} Returns the square length of the vector
   */
  lengthSq() {
    return this.dot(this);
  }

  /**
   * Calculates the length of the vector
   * @return {number} Returns the length of the vector
   */
  length() {
    return Math.sqrt(this.lengthSq());
  }

  /**
   * Returns the magnitude of the passed vector.
   * Sort of like the vector's speed.
   * A vector with a larger x or y will have a larger magnitude.
   * @return {number}
   */
  magnitude() {
    return this.length();
  }

  /**
   * Sets this vector components to the maximum value when
   * compared to the passed vector's components
   * @param {Vec} v The vector to compare to
   * @return {Vec} Returns itself.
   */
  max(v) {
    if (this.x < v.x) {
      this.x = v.x;
    }

    if (this.y < v.y) {
      this.y = v.y;
    }

    return this;
  }

  /**
   * Sets this vector components to the minimum value when
   * compared to the passed vector's components
   * @param {Vec} v The vector to compare to
   * @return {Vec} Returns itself.
   */
  min(v) {
    if (this.x > v.x) {
      this.x = v.x;
    }

    if (this.y > v.y) {
      this.y = v.y;
    }

    return this;
  }

  /**
   * Multiplies the x and y components of this vector by a scalar value
   * @param {number} s The value to multiply by
   * @return {Vec} Returns itself.
   */
  multiplyBy(s) {
    this.x *= s;
    this.y *= s;

    return this;
  }

  /**
   * Negates this vector (multiplies by -1)
   * @return {Vec} Returns itself.
   */
  negate() {
    return this.multiplyBy(-1);
  }

  /**
   * Normalizes this vector (divides by its length)
   * @return {Vec} Returns the normalized vector
   */
  normalize() {
    return this.divideBy(this.length());
  }

  /**
   *
   * @param {Vec} v
   */
  plusEq(v) {
    this.x += v.x;
    this.y += v.y;
  }

  /**
   * Project this vector on to another vector.
   * @param {Vec} v The vector to project onto.
   * @return {Vec} Returns itself.
   */
  project(v) {
    let amt = this.dot(v) / v.lengthSq();
    this.x = amt * v.x;
    this.y = amt * v.y;

    return this;
  }

  /**
   * Project this vector onto a vector of unit length.
   * @param {Vec} v The unit vector to project onto.
   * @return {Vec} Returns itself.
   */
  projectN(v) {
    let amt = this.dot(v);
    this.x = amt * v.x;
    this.y = amt * v.y;

    return this;
  }

  /**
   * Reflect this vector on an arbitrary axis.
   * @param {Vec} axis The vector representing the axis.
   * @return {Vec} Returns itself.
   */
  reflect(axis) {
    let x = this.x,
      y = this.y;
    this.project(axis).multiplyBy(2);
    this.x -= x;
    this.y -= y;

    return this;
  }

  /**
   * Reflect this vector on an arbitrary axis (represented by a unit vector)
   * @param {Vec} axis The unit vector representing the axis.
   * @return {Vec} Returns itself.
   */
  reflectN(axis) {
    let x = this.x,
      y = this.y;
    this.projectN(axis).multiplyBy(2);
    this.x -= x;
    this.y -= y;

    return this;
  }

  /**
   * Rotates the vector by an arbitrary angle
   * around an arbitrary point in space
   * @param {number} angle The angle in radians to rotate by
   * @param {Vec} anchor The anchor point to rotate around
   * @return {Vec} Returns itself.
   */
  rotateAround(angle, anchor) {
    let dist = anchor.distanceTo(this);
    return this.set(
      anchor.x + (dist * Math.cos(angle)),
      anchor.y + (dist * Math.sin(angle))
    );
  }

  /**
   * Rotate the Vec clockwise
   * @param {number} angle The angle in radians to rotate by
   * @return {Vec}
   */
  rotate(angle) {
    let X = this.x * Math.cos(angle) + this.y * Math.sin(angle),
      Y = -this.x * Math.sin(angle) + this.y * Math.cos(angle);

    return new Vec(X, Y);
  }

  /**
   * Round the Vector
   * @return {Vec}
   */
  round() {
    this.x = Math.round(this.x);
    this.y = Math.round(this.y);

    return this;
  }

  /**
   * Scale the Vector
   * @param {number} scale
   * @return {Vec}
   */
  scale(scale) {
    this.x *= scale;
    this.y *= scale;

    return this;
  }

  /**
   * Set the Vec properties
   * @param {number} x
   * @param {number} y
   * @param {number} vx
   * @param {number} vy
   * @param {number} ax
   * @param {number} ay
   */
  set(x, y, vx, vy, ax, ay) {
    this.x = x;
    this.y = y;
    this.vx = vx || this.vx;
    this.vy = vy || this.vy;
    this.ax = ax || this.ax;
    this.ay = ay || this.ay;

    return this;
  }

  /**
   * Sets the length of the vector
   * @param {number} l The length to set this vector to
   * @return {Vec}
   */
  setLength(l) {
    let oldLength = this.length();

    if (oldLength !== 0 && l !== oldLength) {
      this.multiplyBy(l / oldLength);
    }

    return this;
  }

  /**
   * Subtracts a vector from this one
   * @param {Vec} v The vector to subtract from this one
   * @return {Vec}
   */
  subByVec(v) {
    if (typeof v === 'undefined') {
      console.log("Can't sub a vector that is not a vector.");
    }
    return new Vec(this.x - v.x, this.y - v.y);
  }

  /**
   * Subtracts two vectors from each other and
   * stores the result in this vector
   * @param {Vec} a
   * @param {Vec} b
   * @return {Vec} Returns itself.
   */
  subByVecs(a, b) {
    this.x = a.x - b.x;
    this.y = a.y - b.y;

    return this;
  }

  /**
   * Returns an array with the components of this vector as the elements
   * @return {Array}
   */
  toArray() {
    return [Math.round(this.x), Math.round(this.y)];
  }

  /**
   * Convert coords to string
   * @return {string}
   */
  toString() {
    return this.toArray().join(',');
  }

  /**
   * Returns the unit vector for `vector`.
   * A unit vector points in the same direction as the original, but has
   * a magnitude of 1.
   * It's like a direction with a speed that is the same as all other
   * unit vectors.
   * @return {Vec}
   */
  unitVector() {
    return this.divideBy(this.length());
  }
}


/**
 *
 * @param {Vec} v1
 * @param {Vec} v2
 * @return {Vec}
 */
Vec.add = function (v1, v2) {
  return new Vec(v1.x + v2.x, v1.y + v2.y);
};

/**
 *
 * @param {Vec} v1
 * @param {Vec} v2
 * @return {number}
 */
Vec.angleBetween = function (v1, v2) {
  var dotValue = v1.dot(v2);

  return Math.acos(dotValue / (v1.mag() * v2.mag()));
};

/**
 *
 * @param {Vec} v1
 * @param {Vec} v2
 * @return {number}
 */
Vec.distance = function (v1, v2) {
  var dx = v1.x - v2.x,
    dy = v1.y - v2.y;

  return Math.sqrt(dx * dx + dy * dy);
};

/**
 *
 * @param {Vec} v1
 * @param {Vec} v2
 * @return {number}
 */
Vec.distanceSq = function (v1, v2) {
  var dx = v1.x - v2.x,
    dy = v1.y - v2.y;

  return dx * dx + dy * dy;
};

/**
 *
 * @param {Vec} v
 * @param {number} value
 * @return {Vec}
 */
Vec.div = function (v, value) {
  return new Vec(v.x / value, v.y / value);
};

/**
 *
 * @param {Vec} v1
 * @param {Vec} v2
 * @return {boolean}
 */
Vec.equal = function (v1, v2) {
  return (v1.x === v2.x && v1.y === v2.y);
};

/**
 *
 * @param {number} angle
 * @param {number} magnitude
 * @return {Vec}
 */
Vec.fromAngle = function (angle, magnitude) {
  if (magnitude === undefined) {
    magnitude = 1;
  }
  var newVector = new Vec();
  newVector.x = Math.cos(angle) * magnitude;
  newVector.y = Math.sin(angle) * magnitude;

  return newVector;
};

/**
 *
 * @param {b2.Vec} b2Vec
 * @param {boolean} multiply30
 * @return {Vec}
 */
Vec.fromb2Vec = function (b2Vec, multiply30) {
  if (multiply30 === undefined) {
    multiply30 = false;
  }
  if (multiply30) {
    return new Vec(b2Vec.x * 30, b2Vec.y * 30);
  } else {
    return new Vec(b2Vec.x, b2Vec.y);
  }
};

/**
 *
 * @param {Vec} point
 * @param {Vec} linePtA
 * @param {Vec} linePtB
 */
Vec.getNormalPoint = function (point, linePtA, linePtB) {
  var pa = Vec.sub(point, linePtA),
    ba = Vec.sub(linePtB, linePtA);
  ba.normalize().mult(pa.dot(ba));

  return Vec.add(linePtA, ba);
};

/**
 * Get a point at a % point between this Vec and another
 * @param {Vec} v1
 * @param {Vec} v2
 * @return {Vec} .
 */
Vec.getVectorBetween = function (v1, v2) {
  let x = v2.x - v1.x,
    y = v2.y - v1.y;

  return new Vec(x, y);
};

/**
 *
 * @param {Vec} v1
 * @param {Vec} v2
 * @param {number} fraction
 * @return {Vec}
 */
Vec.lerp = function (v1, v2, fraction) {
  var tempV = Vec.sub(v2, v1);
  tempV.mult(fraction);
  tempV.addVecTo(v1);

  return tempV;
};

/**
 *
 * @param {Vec} v1
 * @param {number} value
 * @return {Vec}
 */
Vec.mult = function (v1, value) {
  return new Vec(v1.x * value, v1.y * value);
};

/**
 *
 * @param {number} xmin
 * @param {number} xmax
 * @param {number} ymin
 * @param {number} ymax
 * @return {Vec}
 */
Vec.random = function (xmin, xmax, ymin, ymax) {
  var result = new Vec();
  result.x = Utility.Maths.map(Math.random(), 0, 1, xmin, xmax);
  result.y = Utility.Maths.map(Math.random(), 0, 1, ymin, ymax);

  return result;
};

/**
 *
 * @param {Vec} v
 * @param {number} angle
 */
Vec.setAngle = function (v, angle) {
  var mag = v.mag();
  v.x = Math.cos(angle) * mag;
  v.y = Math.sin(angle) * mag;
};

/**
 *
 * @param {Vec} v1
 * @param {Vec} v2
 * @return {Vec}
 */
Vec.sub = function (v1, v2) {
  return new Vec(v1.x - v2.x, v1.y - v2.y);
};

/**
 * Get a Vec between this Vec and another
 * @param {Vec} v1
 * @param {Vec} v2
 * @return {Vec} .
 */
Vec.vectorBetween = function (v1, v2) {
  let x = v2.x - v1.x,
    y = v2.y - v1.y;

  return new Vec(x, y);
};
