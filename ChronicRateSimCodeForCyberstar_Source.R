#-- functions to source for simulations that estimate carriage rate --#
#
#static.graph <- graph.list[[1]]
#prop.chronic = .2
#P.mort.given.chronic.contact = .8 
#max.ts = 8
#edgelists.nozeros = edge.list[[1]]
#
#-- 1) function to transmit disease through graph --#
disease.through.graph <- function(static.graph, prop.chronic,
																	P.mort.given.chronic.contact, max.ts,
																	edgelists.nozeros){
#
#static.graph <- graph.list[[1]]
#prop.chronic <- .2
#P.mort.given.chronic.contact <- .8
#max.ts <- 8
#edgelists.nozeros <- edgelist.nozeros[[1]]
#
#

still.alive <- edgelist.ofinterest <- static.graph.list <- el.ofinterest <- chronic.nodes  <- vector("list", length = max.ts)
edgelist.ofinterest[[1]] <- edgelists.nozeros
el.ofinterest[[1]] <- cbind(edgelist.ofinterest[[1]]$Ind1,
														edgelist.ofinterest[[1]]$Ind2)

static.graph.list[[1]] <- static.graph
graph.edgelist(el.ofinterest[[1]], directed = F)
prop.chronic <- prop.chronic 
P.mort.given.chronic.contact <-		 P.mort.given.chronic.contact 
surviving <- rep(NA, max.ts)

#-- main loop -#
for(j in 1:max.ts){
if(j == 1){
#-- in first timestep, select initially infected nodes --#
chronic.nodes[[j]] <- rbinom(length(V(static.graph.list[[1]])$name), 1,
														 prop.chronic)
chronic.names <- V(static.graph.list[[1]])$name[chronic.nodes[[j]] == 1]
new.I.status <- rep(NA, length(V(static.graph.list[[j]])$name))
prob.list <- vector("list", length(levels(factor(V(static.graph.list[[j]])$name))))
for(i in 1:length(levels(factor(V(static.graph.list[[j]])$name)))){
#-- check to			see if this node is already infected --#
#-- if this node ISN'T alre			ady infected, do this: --#
if(! (V(static.graph.list[[j]])$name[i] %in% chronic.names)){ 
#-- extract all edges linking to node i --#
					k <- subset(edgelist.ofinterest[[j]], as.character(Ind1) == V(static.graph)$name[i] | as.character(Ind2) == V(static.graph)$name[i])
#-- label all nodes with their infection status --#
					m <- subset(k, as.character(Ind1) %in% chronic.names | as.character(Ind2)%in% chronic.names)
#-- sum edgeweights to infected nodes													, multiply by
#-- P.mort.given.chronic.contact --#
		prob.become.infected <- 1 - exp(-sum(m$edgeweights) *
																	P.mort.given.chronic.contact)
		prob.list[[i]] <- prob.become.infected
		infection.trial <- rbinom(1, 1, prob.become.infected)
		new.I.status[i] <- ifelse(infection.trial == 1, "I", "S")
			} else {
				new.I.status[i] <- "Dead"
			}
		}
		still.alive[[j]] <- subset(data.frame(cbind(V(static.graph)$name, new.I.status)), new.I.status != "Dead")
		names(still.alive[[j]]) <- c("Ind", "Status")
		#-- end if(j = 1)		== T --# 
		} else { #-- j > 1 == T --#
		prev.lengths <- rep(NA, j - 1)
			for(b in 1:(j - 1)){
					prev.lengths[b] <- length(V(static.graph.list[[b]])$name)
			}
			if(! is.na(table(prev.lengths == 0)["TRUE"])){
													surviving[j] <- 0
			} else {
			#-- in later timesteps, chronic.nodes are nodes currently in I --#
				new.I.status <- rep(NA, length(V(static.graph.list[[j]])$name))
				chronic.nodes[[j]] <- subset(still.alive[[j - 1]],
													 as.character(still.alive[[j - 1]]$Status) == "I")$Ind 
					for(i in 1:length(levels(factor(V(static.graph.list[[j]])$name)))){
													#-- check to see if this node is already infected --#
													#-- if this node ISN'T already infected, do this: --#
		if(!(levels(factor(V(static.graph.list[[j]])$name))[i] %in% chronic.names)){
													#-- extract all edges linking to node i --#
				k <- subset(edgelist.ofinterest[[j]], as.character(Ind1) == V(static.graph)$name[i] | as.character(Ind2) == V(static.graph)$name[i])
							#-- label all nodes with their infection status --#
				m <- subset(k, as.character(Ind1) %in% chronic.names | as.character(Ind2)
				%in% chronic.names)
							#-- sum edgeweights to infected nodes, multiply by
							#-- P.mort.given.chronic.contact --#
				prob.become.infected <- 1 - exp(-sum(m$edgeweights) *
																	P.mort.given.chronic.contact)
				infection.trial <- rbinom(1, 1, prob.become.infected)
				new.I.status[i] <- ifelse(infection.trial == 1, "I", "S")
				} else {
						new.I.status[i] <- "Dead"
				}
		}
		still.alive[[j]] <- subset(data.frame(cbind(V(static.graph.list[[j]])$name, new.I.status)), new.I.status != "Dead")
#		still.alive[[j]] <- subset(data.frame(cbind(levels(factor(edgelist.ofinterest[[j]]$Ind1)), new.I.status)), new.I.status != "Dead")
		names(still.alive[[j]]) <- c("Ind", "Status")
			}
		}
		#-- update graph by removing all Dead nodes --#
		edgelist.ofinterest[[j + 1]] <- subset(edgelist.ofinterest[[j]], as.character(Ind1) %in% as.character(still.alive[[j]]$Ind) & as.character(Ind2) %in% as.character(still.alive[[j]]$Ind))
		el <- cbind(as.character(edgelist.ofinterest[[j + 1]]$Ind1),
					as.character(edgelist.ofinterest[[j + 1]]$Ind2))
		static.graph.list[[j + 1]] <- graph.edgelist(el, directed = F)
		E(static.graph.list[[j + 1]])$weight <- edgelist.ofinterest[[j +
																																 1]]$edgeweights
#}
	surviving[j] <- length(levels(factor(V(static.graph.list[[j]])$name)))
	}
	return(surviving)
}
#
##-- specify path to source file --#
#root.path <- "~/work/CompoChronicSims/" 
#sourcepath <-
#	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Code/SimsToEstChronicRate/"
#write.path <- "~/work/CompoChronicSims/Output/"
#indices <- 1:8
#
#	#-- specify path to data files --#
#graphpath <-
#	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/GraphObjectData_27Sept2013/StaticGraphList" 
###
#edgelist.path <-
#	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/GraphObjectData_27Sept2013/EdgelistsNozeros"
###
##-- read in source and data files --#
##source(paste(rootpath, "Code/ChronicRateSimCodeForCyberstar_Source.R", sep = ""))
##source(paste(sourcepath, "ChronicRateSimCodeForCyberstar_Source.R", sep = ""))
##source(paste(sourcepath, "ChronicRateSimCodeForCyberstar_Source_V2.R", sep = ""))
##graph.list.full <- dget(paste(root.path, "Data/StaticGraphList", sep = ""))
#graph.list.full <- dget(paste(graphpath, sep = ""))
#edgelist.nozeros.full <- dget(paste(edgelist.path, sep = ""))
##edgelist.nozeros.full <- dget(paste(root.path, "Data/EdgelistsNoZeros"))
##graph.list.full <- dget(paste(graphpath, sep = ""))
##-- subset down to smaller group of elements for batching to Cyberstar --#
#graph.list <- graph.list.full[indices]
#edgelist.nozeros <- edgelist.nozeros.full[indices]
#
#-- call required packages --#
#require(igraph)
#
#test1 <- disease.through.graph(graph.list[[1]], prop.chronic = .2,
#															 P.mort.given.chronic.contact = .8, max.ts = 8,
#															 edgelists.nozeros = edgelist.nozeros[[1]])
#
#require(igraph)
#bottompath <-
#	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/GraphObjectData_27Sept2013/"
#
#graph.list <- dget(paste(bottompath, "StaticGraphList", sep = ""))
#edge.list <- dget(paste(bottompath, "EdgelistsNozeros", sep = ""))
#
#test1 <- disease.through.graph(graph.list[[1]], prop.chronic = .2,
#															 P.mort.given.chronic.contact = .8, max.ts = 8,
#															 edgelists.nozeros = edge.list[[1]])

#-- batch function --#
#-- replicate disease.through.graph many times. --#
batch.fun <- function(static.graph = static.graph, prop.chronic,
		P.mort.given.chronic.contact, max.ts,
		edgelists.nozeros, reps){
		lamb.mort.vec <- test <- tot.ewe.vec <- rep(NA, reps)
			for(r in 1:reps){
			#lamb.mort.vec[r] <- disease.through.graph(static.graph = static.graph, 
			#prop.chronic = prop.chronic,
			#P.mort.given.chronic.contact
			#= P.mort.given.chronic.contact, 
			#max.ts = max.ts, 
			#edgelists.nozeros =
			#edgelists.nozeros)[max.ts]
				test[r] <- try(disease.through.graph(static.graph = static.graph, 
											prop.chronic = prop.chronic,
											P.mort.given.chronic.contact = P.mort.given.chronic.contact, 
											max.ts = max.ts, 
											edgelists.nozeros = edgelists.nozeros)[max.ts])
				lamb.mort.vec[r] <- ifelse(class(test[r]) != "try-error", test[r], NA)
				tot.ewe.vec[r] <- length(levels(factor(V(static.graph)$name))) 
			print(r)
		}
		out.mat <- cbind(lamb.mort.vec, tot.ewe.vec)
		return(out.mat)
}

param.space.fun <- function(min.Pmort, max.Pmort, mortstep, min.Pchron,
														max.Pchron, chronstep, static.graph, reps,
														edgelists.nozeros, max.ts = max.ts){
	mort.steps <- seq(min.Pmort, max.Pmort, by = mortstep)
	chron.steps <- seq(min.Pchron, max.Pchron, by = chronstep)
	param.mat <- expand.grid(mort.steps, chron.steps)
	out.list <- vector("list", dim(param.mat)[1])
	for(p in 1:dim(param.mat)[1]){
		out.list[[p]] <- batch.fun(static.graph, prop.chronic = param.mat[p, 1],
															 P.mort.given.chronic.contact = param.mat[p, 2],
															 max.ts = max.ts, edgelists.nozeros =
															 edgelists.nozeros, reps = reps)
	}
	return.list <- list(param.mat = param.mat, out.list = out.list)
	return(return.list)
}





