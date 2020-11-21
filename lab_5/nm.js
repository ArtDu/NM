function solve(a, b, c, alpha, beta, gamma, delta, x, h, t, thau, scheme, approx) {
    function psi(x) {
        return Math.sin(x);
    }

    function f(x, t) {
        return 0;
    }


    function phi0(t) {
        return Math.exp((c - a) * t);
    }

    function phi1(t) {
        return Math.exp((c - a) * t);
    }


    let sigma = a ** 2 * thau / h ** 2;
    document.getElementById('sigma').textContent = "σ = " + sigma.toString();
    if (scheme === 1) {
        console.log('sigma must be <= 0.5', sigma);
        if (sigma >= 0.5) {
            alert('Явная схема неустойчива при заданных параметрах. Решения будут неточными. Следует, например, увеличить разбиение по t.')
        }
    }

    let n = x.length - 1;
    let U = zeroes(t.length, x.length);

    for (let j = 0; j < x.length; j++)
        U[0][j] = psi(x[j]);

    if (scheme === 1) {
        // Явная схема
        for (let k = 1; k < t.length; k++) {
            for (let j = 1; j < x.length - 1; j++) {
                U[k][j] = (U[k - 1][j + 1] * (a * thau / h ** 2 + b * thau / 2 / h)
                    + U[k - 1][j] * (-2 * a * thau / h ** 2 + c * thau + 1)
                    + U[k - 1][j - 1] * (a * thau / h ** 2 - b * thau / 2 / h)
                    + thau * f(x[j], t[k]));
            }
            if (approx === 1) {
                // u[k + 1][0] = (h * phi0(t[k + 1], a, b, c) - alpha * u[k + 1][1]) / (h * beta - alpha);
                U[k][0] = (phi0(t[k]) - alpha / h * U[k][1]) / (beta - alpha / h);
                U[k][n] = (phi1(t[k]) + gamma / h * U[k][n - 1]) / (delta + gamma / h);
                // u[k + 1][n] = (h * phi1(t[k + 1], a, b, c) + gamma * u[k + 1][n - 1]) / (h * delta + gamma);
            } else if (approx === 2) {
                U[k][0] = (((2 * alpha * a / h / (2 * a - h * b)) * U[k][1] + (alpha * h / thau / (2 * a - h * b)) * U[k - 1][0]
                    + (alpha * h / (2 * a - h * b)) * f(0, t[k]) - phi0(t[k]))
                    / ((2 * alpha * a / h / (2 * a - h * b)) + (alpha * h / thau / (2 * a - h * b)) - (alpha * h / (2 * a - h * b)) * c - beta));
                U[k][n] = (((2 * gamma * a / h / (2 * a + h * b)) * U[k][n - 1] + (gamma * h / thau / (2 * a + h * b)) * U[k - 1][n] + (gamma * h * c / (2 * a + h * b)) * f(l, t[k]) + phi1(t[k]))
                    / ((2 * gamma * a / h / (2 * a + h * b)) + (gamma * h / thau / (2 * a + h * b)) - (gamma * h * c / (2 * a + h * b)) * c + delta));
            } else if (approx === 3) {
                U[k][0] = ((phi0(t[k]) - U[k][1] * 2 * alpha / h + U[k][2] * alpha / 2 / h)
                    / (beta - 3 * alpha / 2 / h));
                U[k][n] = ((phi1(t[k]) + U[k][n - 1] * 2 * gamma / h - U[k][n - 2] * gamma / 2 / h)
                    / (delta + 3 * gamma / 2 / h))
            }
        }
        return U;
    }
    if (scheme === 2) {
        // Неявная схема
        for (let k = 1; k < t.length; k++) {

            let diag_a = zeroes(x.length);
            let diag_b = zeroes(x.length);
            let diag_c = zeroes(x.length);
            let column = zeroes(x.length);

            if (approx === 1) {
                diag_a[0] = 0;
                diag_b[0] = -alpha / h + beta;
                diag_c[0] = alpha / h;
                column[0] = phi0(t[k]);

                diag_a[n] = -gamma / h;
                diag_b[n] = gamma / h + delta;
                diag_c[n] = 0;
                column[n] = phi1(t[k]);
            }
            if (approx === 2) {
                let d = 2 * a - h * b;
                let K1 = 2 * alpha * a / h / d;
                let K2 = alpha * h / thau / d;
                let K3 = alpha * h / d;
                diag_a[0] = 0;
                diag_b[0] = K1 + K2 - K3 * c - beta;
                diag_c[0] = -K1;
                column[0] = K2 * U[k - 1][0] + K3 * f(0, t[k]) - phi0(t[k]);

                d = 2 * a + h * b;
                K1 = 2 * gamma * a / h / d;
                K2 = gamma * h / thau / d;
                K3 = gamma * h / d;
                diag_a[n] = -K1;
                diag_b[n] = (K1 + K2 - K3 * c + delta);
                diag_c[n] = 0;
                column[n] = K2 * U[k - 1][n] + K3 * f(x[n], t[k]) + phi1(t[k]);
            }
            if (approx === 3) {
                diag_a[0] = -3 * alpha / 2 / h + beta;
                diag_b[0] = 2 * alpha / h;
                diag_c[0] = -alpha / 2 / h;
                column[0] = phi0(t[k]);

                diag_a[n] = gamma / 2 / h;
                diag_b[n] = -2 * gamma / h;
                diag_c[n] = 3 * gamma / 2 / h + delta;
                column[n] = phi1(t[k]);
            }


            for (let j = 1; j < x.length - 1; ++j) {
                diag_a[j] = -thau * a / h ** 2 + thau * b / 2 / h;
                diag_b[j] = 1 + 2 * thau * a / h ** 2 - thau * c;
                diag_c[j] = -thau * a / h ** 2 - thau * b / 2 / h;
                column[j] = U[k - 1][j] + thau * f(x[j], t[k]);
            }

            if (approx === 3) {
                let const_k = diag_c[0] / diag_c[1];
                diag_a[0] -= const_k * diag_a[1];
                diag_b[0] -= const_k * diag_b[1];
                column[0] -= const_k * column[1];
                diag_c[0] = diag_b[0];
                diag_b[0] = diag_a[0];
                diag_a[0] = 0;

                const_k = diag_a[n] / diag_a[n - 1];
                diag_b[n] -= const_k * diag_b[n - 1];
                diag_c[n] -= const_k * diag_c[n - 1];
                column[n] -= const_k * column[n - 1];
                diag_a[n] = diag_b[n];
                diag_b[n] = diag_c[n];
                diag_c[n] = 0;
            }

            let xx = tridiagonal(diag_a, diag_b, diag_c, column);
            for (let j = 0; j <= x.length; ++j) {
                U[k][j] = xx[j];
            }
        }
        return U;
    }
    if (scheme === 3) {
        for (let k = 1; k < t.length; k++) {

            let diag_a = zeroes(x.length);
            let diag_b = zeroes(x.length);
            let diag_c = zeroes(x.length);
            let column = zeroes(x.length);

            if (approx === 1) {
                diag_a[0] = 0;
                diag_b[0] = -alpha / h + beta;
                diag_c[0] = alpha / h;
                column[0] = phi0(t[k]);

                diag_a[n] = -gamma / h;
                diag_b[n] = gamma / h + delta;
                diag_c[n] = 0;
                column[n] = phi1(t[k]);
            }
            if (approx === 2) {
                let d = 2 * a - h * b;
                let coeff1 = 2 * alpha * a / h / d;
                let coeff2 = alpha * h / thau / d;
                let coeff3 = alpha * h / d;
                diag_a[0] = 0;
                diag_b[0] = coeff1 + coeff2 - coeff3 * c - beta;
                diag_c[0] = -coeff1;
                column[0] = coeff2 * U[k - 1][0] + coeff3 * f(0, t[k]) - phi0(t[k]);

                d = 2 * a + h * b;
                coeff1 = 2 * gamma * a / h / d;
                coeff2 = gamma * h / thau / d;
                coeff3 = gamma * h / d;
                diag_a[n] = -coeff1;
                diag_b[n] = (coeff1 + coeff2 - coeff3 * c + delta);
                diag_c[n] = 0;
                column[n] = coeff2 * U[k - 1][n] + coeff3 * f(x[n], t[k]) + phi1(t[k]);
            }
            if (approx === 3) {
                diag_a[0] = -3 * alpha / 2 / h + beta;
                diag_b[0] = 2 * alpha / h;
                diag_c[0] = -alpha / 2 / h;
                column[0] = phi0(t[k]);

                diag_a[n] = gamma / 2 / h;
                diag_b[n] = -2 * gamma / h;
                diag_c[n] = 3 * gamma / 2 / h + delta;
                column[n] = phi1(t[k]);
            }


            for (let j = 1; j < x.length - 1; ++j) {
                diag_a[j] = (-thau * a / h ** 2 + thau * b / 2 / h);
                diag_b[j] = (2 + 2 * thau * a / h ** 2 - thau * c);
                diag_c[j] = (-thau * a / h ** 2 - thau * b / 2 / h);
                column[j] = (U[k - 1][j - 1] * (thau * a / h ** 2 - thau * b / 2 / h)
                    + U[k - 1][j] * (2 - 2 * thau * a / h ** 2 + thau * c)
                    + U[k - 1][j + 1] * (thau * a / h ** 2 + thau * b / 2 / h)
                    + thau * (f(x[j], t[k]) + f(x[j], t[k - 1])));
            }

            if (approx === 3) {
                let const_k = diag_c[0] / diag_c[1];
                diag_a[0] -= const_k * diag_a[1];
                diag_b[0] -= const_k * diag_b[1];
                column[0] -= const_k * column[1];
                diag_c[0] = diag_b[0];
                diag_b[0] = diag_a[0];
                diag_a[0] = 0;

                const_k = diag_a[n] / diag_a[n - 1];
                diag_b[n] -= const_k * diag_b[n - 1];
                diag_c[n] -= const_k * diag_c[n - 1];
                column[n] -= const_k * column[n - 1];
                diag_a[n] = diag_b[n];
                diag_b[n] = diag_c[n];
                diag_c[n] = 0;

            }

            let xx = tridiagonal(diag_a, diag_b, diag_c, column);
            for (let j = 0; j <= x.length; ++j) {
                U[k][j] = xx[j];
            }
        }
        return U;
    }
}


function zeroes(rows, cols = 0) {
    let i;
    let arr = [];

    if (cols === 0) {
        for (i = 0; i < rows; i++)
            arr.push(0);
    } else {
        for (i = 0; i < rows; i++)
            arr.push([]);

        for (i = 0; i < rows; i++)
            for (let j = 0; j < cols; j++)
                arr[i].push(0);
    }
    return arr;
}

function arange(from, to, step) {
    let arr = [];
    while (from < to) {
        arr.push(from);
        from += step;
    }
    return arr;
}



function tridiagonal(a, b, c, d) {
    let n = d.length;
    let p = [-c[0] / b[0]];
    let q = [d[0] / b[0]];

    for (let i = 1; i < n; i++) {
        p.push(-c[i] / (b[i] + a[i] * p[i - 1]));
        q.push((d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]));
    }
    let x = zeroes(n);
    x[n - 1] = q[n - 1];
    for (let i = n - 2; i >= 0; --i) {
        x[i] = p[i] * x[i + 1] + q[i];
    }

    return x
}

function exact(x, t) {
    return Math.exp((c - a) * t) * Math.sin(x);

}

function analytic_solve(x_list, t, a, b, c) {
    let solved = [];
    x_list.forEach(item =>
        solved.push(exact(item, t, a, b, c))
    );
    return solved;
}

function get_diff(x, y) {
    let diff = [];
    for (let i = 0; i < Math.min(x.length, y.length); i++) {
        diff.push(Math.abs(x[i] - y[i]));
    }
    return diff;
}
