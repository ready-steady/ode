// Package rk4 provides an integrator of system of ordinary differential
// equations based on the fourth-order Runge–Kutta method.
//
// https://en.wikipedia.org/wiki/Runge–Kutta_methods
package rk4

// Integrator is an integrator.
type Integrator struct {
	config Config
}

// New creates a new integrator.
func New(config *Config) (*Integrator, error) {
	if err := config.verify(); err != nil {
		return nil, err
	}
	return &Integrator{config: *config}, nil
}

// Compute integrates the system of differential equations dy/dx = f(x, y). See
// Integrator.Compute in the parent package.
//
// The solution is returned at a number of equidistant points starting from and
// including x0 = xs[0]. The final point is the closest point to the last
// element of xs with respect to the integration step.
func (self *Integrator) Compute(dydx func(float64, []float64, []float64),
	y0 []float64, xs []float64) ([]float64, []float64, error) {

	nd, nx := len(y0), len(xs)

	z := make([]float64, nd)

	f := make([]float64, 4*nd)
	f1 := f[0*nd : 1*nd]
	f2 := f[1*nd : 2*nd]
	f3 := f[2*nd : 3*nd]
	f4 := f[3*nd : 4*nd]

	x0, xend := xs[0], xs[nx-1]

	h := self.config.Step

	ns := int((xend-x0)/h+0.5) + 1

	// Done with the first point.
	ys := make([]float64, ns*nd)
	copy(ys, y0)

	for k, x, y := 1, x0, y0; k < ns; k++ {
		// Step 1
		dydx(x, y, f1)

		// Step 2
		for i := 0; i < nd; i++ {
			z[i] = y[i] + h*f1[i]/2
		}
		dydx(x+h/2, z, f2)

		// Step 3
		for i := 0; i < nd; i++ {
			z[i] = y[i] + h*f2[i]/2
		}
		dydx(x+h/2, z, f3)

		// Step 4
		for i := 0; i < nd; i++ {
			z[i] = y[i] + h*f3[i]
		}
		dydx(x+h, z, f4)

		ynew := ys[k*nd:]
		for i := 0; i < nd; i++ {
			ynew[i] = y[i] + h*(f1[i]+2*f2[i]+2*f3[i]+f4[i])/6
		}

		x += h
		y = ynew
	}

	return ys, xs, nil
}
