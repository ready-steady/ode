package dopri

import (
	"errors"
)

// Config is the configuration of an integrator.
type Config struct {
	// The initial step of integration.
	TryStep float64
	// The maximal step of integration.
	MaxStep float64
	// The absolute error tolerance.
	AbsError float64
	// The relative error tolerance.
	RelError float64
}

// DefaultConfig returns the default configuration of an integrator.
func DefaultConfig() *Config {
	return &Config{
		TryStep:  0,
		MaxStep:  0,
		AbsError: 1e-6,
		RelError: 1e-3,
	}
}

func (c *Config) verify() error {
	if c.TryStep < 0 {
		return errors.New("the initial step should be nonnegative")
	}
	if c.MaxStep < 0 {
		return errors.New("the maximal step should be nonnegative")
	}
	if c.AbsError <= 0 {
		return errors.New("the absolute error tolerance should be positive")
	}
	if c.RelError <= 0 {
		return errors.New("the relative error tolerance should be positive")
	}

	return nil
}
