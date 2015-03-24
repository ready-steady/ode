// Package ode provides algorithms for integrating systems of ordinary
// differential equations.
package ode

type Integrator interface {
	// Compute integrates the system of differential equations dy/dx = f(x, y).
	//
	// The input function dydx(x, y, f) evaluates f(x, y) for a given x and y in
	// its first and second arguments and stores the result in its third
	// argument. The interval of integration is [x0, xend] where x0 and xend are
	// the first and last entries of xs, respectively. The initial condition is
	// y0, which correspond to x0.
	Compute(dydx func(float64, []float64, []float64), y0 []float64,
		xs []float64) ([]float64, []float64, error)
}
