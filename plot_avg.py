# Trying to explain why sample 350 is magic. Does it fall near to a clock rising edge?
# Answer is "NO", since there is not clock period in the traces... maybe 62 samples, but irregular (see below), sometimes 47/48 as well.

import numpy as np

D = 700
mean_trace = np.zeros(D, dtype=float)
nb_traces = 0

with open( 'Attack_traces_traces.txt', 'r' ) as f:
	for trace in f:
		v = trace.split(',')
		assert len(v) == D
		mean_trace += np.array( [ float(c) for c in v ] ) # map( float, v )
		nb_traces += 1

with open( 'avg.txt', 'w' ) as f:
	for sample in mean_trace / nb_traces:
		f.write( "%s\n"%sample )

# Maxima, and delta:
#
#  27
#  90	63
# 137	47
# 200	63
# 262	62
# 325	63
# 388	63
# 436	48
# 498	62
# 561	63
# 623	62
