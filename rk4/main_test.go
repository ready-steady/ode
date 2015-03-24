package rk4

import (
	"testing"

	"github.com/ready-steady/support/assert"
)

func TestComputeSimple(t *testing.T) {
	fixture := &fixtureSimple

	integrator, _ := New(fixture.configure())
	ys, _, _ := integrator.Compute(fixture.dydx, fixture.y0, fixture.xs)

	assert.Equal(ys, fixture.ys, t)
}

func TestComputeKraichnanOrszag(t *testing.T) {
	fixture := &fixtureKraichnanOrszag

	integrator, _ := New(fixture.configure())
	ys, _, _ := integrator.Compute(fixture.dydx, fixture.y0, fixture.xs)

	assert.EqualWithin(ys, fixture.ys, 5e-14, t)
}
