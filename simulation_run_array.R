for (tree in paste0("tree", 1:4)) {
	for (fc in c(paste0("linear", 0:2), "poly2")) {
		for (rho in c(0.2, 0.5, 0.8)) {
			for (delta in c(1, 3)) {
				for (n in c(100, 200, 500)) {
					for (k in 1:1000) {
						source("simulation_correct.R")
					}
				}
			}
		}
	}
}
