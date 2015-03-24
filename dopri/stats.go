package dopri

// Stats contains information about the work done by an integrator.
type Stats struct {
	Evaluations uint // The number of invocations of the derivative function.
	Rejections  uint // The number of rejected iterations of the algorithm.
	Steps       uint // The number of steps the algorithm has taken.
}
