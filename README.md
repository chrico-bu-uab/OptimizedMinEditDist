# OptimizedMinEditDist
Optimized Minimum Edit Distance (Weighted)

This Jupyter notebook contains code for applying Simple Good-Turing smoothing to spelling correction confusion matrices, and then using those confusion matrices to calculate weighted Damerau-Levenshtein distance. https://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance
This calculation is then optimized to minimize the distance between misspellings and their correct spelling obtained from an augmented version of the following list: https://www.lexico.com/grammar/common-misspellings

See https://pdfs.semanticscholar.org/3c0f/046634f8102c2acb495aaf7f14924c2d4ee7.pdf for a detailed explanation of the smoothing procedure.

The confusion matrix files were borrowed from https://github.com/jbhoosreddy/spellcorrect/ and originally derived in the following paper: https://www.cs.ubc.ca/~carenini/TEACHING/CPSC503-04/spelling90.pdf

The lines on the graphs were created based on the following StackOverflow answer: https://stackoverflow.com/a/43811762/5295786
