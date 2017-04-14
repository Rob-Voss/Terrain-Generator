// Ported from Stefan Gustavson's java implementation
// http://staffwww.itn.liu.se/~stegu/simplexnoise/simplexnoise.pdf
// Sean McCullough banksean@gmail.com

class SimplexNoise {
  constructor(gen) {
    this.rand = gen;
    this.grad3 = [
      [1, 1, 0], [-1, 1, 0], [1, -1, 0], [-1, -1, 0],
      [1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1],
      [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1]
    ];

    this.simplex = [
      [0, 1, 2, 3], [0, 1, 3, 2], [0, 0, 0, 0], [0, 2, 3, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [1, 2, 3, 0],
      [0, 2, 1, 3], [0, 0, 0, 0], [0, 3, 1, 2], [0, 3, 2, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [1, 3, 2, 0],
      [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
      [1, 2, 0, 3], [0, 0, 0, 0], [1, 3, 0, 2], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [2, 3, 0, 1], [2, 3, 1, 0],
      [1, 0, 2, 3], [1, 0, 3, 2], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [2, 0, 3, 1], [0, 0, 0, 0], [2, 1, 3, 0],
      [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
      [2, 0, 1, 3], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [3, 0, 1, 2], [3, 0, 2, 1], [0, 0, 0, 0], [3, 1, 2, 0],
      [2, 1, 0, 3], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [3, 1, 0, 2], [0, 0, 0, 0], [3, 2, 0, 1], [3, 2, 1, 0]
    ];
  }

  setSeed(seed) {
    this.p = [];
    this.rand.seed = seed;

    for (let i = 0; i < 256; i++) {
      this.p[i] = Math.floor(this.rand.nextRange(0, 255));
    }

    this.perm = [];
    for (let i = 0; i < 512; i++) {
      this.perm[i] = this.p[i & 255];
    }
  }

  dot(g, x, y) {
    return g[0] * x + g[1] * y;
  }

  noise(xin, yin) {
    let n0, n1, n2,
      F2 = 0.5 * (Math.sqrt(3.0) - 1.0),
      s = (xin + yin) * F2,
      i = Math.floor(xin + s),
      j = Math.floor(yin + s),
      G2 = (3.0 - Math.sqrt(3.0)) / 6.0,
      t = (i + j) * G2,
      X0 = i - t,
      Y0 = j - t,
      x0 = xin - X0,
      y0 = yin - Y0,
      i1, j1;
    if (x0 > y0) {
      i1 = 1;
      j1 = 0;
    } else {
      i1 = 0;
      j1 = 1;
    }

    let x1 = x0 - i1 + G2,
      y1 = y0 - j1 + G2,
      x2 = x0 - 1.0 + 2.0 * G2,
      y2 = y0 - 1.0 + 2.0 * G2,
      ii = i & 255,
      jj = j & 255,
      gi0 = this.perm[ii + this.perm[jj]] % 12,
      gi1 = this.perm[ii + i1 + this.perm[jj + j1]] % 12,
      gi2 = this.perm[ii + 1 + this.perm[jj + 1]] % 12,
      t0 = 0.5 - x0 * x0 - y0 * y0;
    if (t0 < 0) {
      n0 = 0.0;
    } else {
      t0 *= t0;
      n0 = t0 * t0 * this.dot(this.grad3[gi0], x0, y0);
    }
    let t1 = 0.5 - x1 * x1 - y1 * y1;
    if (t1 < 0) {
      n1 = 0.0;
    } else {
      t1 *= t1;
      n1 = t1 * t1 * this.dot(this.grad3[gi1], x1, y1);
    }
    let t2 = 0.5 - x2 * x2 - y2 * y2;
    if (t2 < 0) {
      n2 = 0.0;
    } else {
      t2 *= t2;
      n2 = t2 * t2 * this.dot(this.grad3[gi2], x2, y2);
    }

    return 70.0 * (n0 + n1 + n2);
  }
}
