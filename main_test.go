package ode

import (
	"testing"

	"github.com/ready-steady/numeric/integration/ode/dopri"
	"github.com/ready-steady/numeric/integration/ode/rk4"
)

func TestIntegrator(t *testing.T) {
	var integrator Integrator

	integrator, _ = dopri.New(dopri.DefaultConfig())
	integrator, _ = rk4.New(&rk4.Config{Step: 42})

	blackbox(integrator)
}

func blackbox(_ interface{}) {}
