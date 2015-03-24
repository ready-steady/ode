package rk4

import (
	"errors"
)

// Config is the configuration of an integrator.
type Config struct {
	// The step of integration.
	Step float64
}

func (c *Config) verify() error {
	if c.Step < 0 {
		return errors.New("the initial step should be nonnegative")
	}

	return nil
}
