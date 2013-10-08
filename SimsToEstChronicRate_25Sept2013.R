#-- Simulation code to estimate chronic carriage rate using empirically --#
#-- structured networks --#

#-- read in assocations data --#
filepath <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/"

data <- read.csv(paste(filepath, "FourPopRelocs_Through2012.csv", sep = ""),
								 header =T)
data <- subset(data, Pop != "" & Pop != "Pop" & Year != "") 

#-- source in Paul's code to build association matrix --#
source.path <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Code/DataCleaning/RelocationsToNetworks/StaticNetworkAssocMat.R"
source(paste(source.path))

#-- subset data into pop-year chunks --#
pops <- levels(factor(data$Pop))
years <- levels(factor(data$Year))
pop.id <- year.id <- rep(NA, length(pops) * length(years))
data.subsets <- list(NA, length(pops) * length(years))
for(i in 1:length(pops)){
	for(j in 1:length(years)){
		data.subsets[[(i-1) * length(years) + j]] <- subset(data, Pop == pops[i] &
																												Year == years[j])
		pop.id[(i-1) * length(years) + j] <- pops[i]
		year.id[(i-1) * length(years) + j] <- years[j]
	}
}

static.graph.list <- edgelists.nozeros <- inds <- edgelists <- output.info <- list(NA, length(pops) *
																											length(years))

require(igraph)
edgeweight.min <- 0
#-- build network --#
#i <- 1
for(i in 1:length(data.subsets)){
if(dim(data.subsets[[i]])[1] == 0 |
	 length(levels(factor(data.subsets[[i]]$EWEID))) == 1){
	edgelists[[i]] <- NA
	inds[[i]] <- NA
	output.info[[i]] <- NA
} else {
	out <- AssocTimeVect(AssocData = data.subsets[[i]])
	inds[[i]] <- out[[1]]
	edgelists[[i]] <- out[[2]]
	edgelists[[i]]$edgeweights <-
		as.numeric(as.character(edgelists[[i]]$TimesTogether)) /
	(as.numeric(as.character(edgelists[[i]]$TotalInd1)) +
	 as.numeric(as.character(edgelists[[i]]$TotalInd2)) -
	 as.numeric(as.character(edgelists[[i]]$TimesTogether)))
	edgelists.nozeros[[i]] <- subset(edgelists[[i]], edgeweights >=
																	 edgeweight.min)
	#-- include edges with weights = 0 --#
	el <- cbind(as.character(edgelists.nozeros[[i]]$Ind1),
							as.character(edgelists.nozeros[[i]]$Ind2))
	#-- generate a list of isolated nodes --#
	connected.nodes <- levels(factor(c(el[, 1], el[, 2])))
	isolated.nodes <- inds[[i]][! (inds[[i]] %in% connected.nodes) == T]
	if(dim(el)[1] == 0){
		static.graph.orig <- graph.edgelist(el, directed = F)
	} else {
	static.graph.orig <- graph.edgelist(el, directed = F)
	static.graph.weighted <- set.edge.attribute(static.graph.orig, "weight",
																							value =
																							edgelists.nozeros[[i]]$edgeweights)
	E(static.graph.weighted)$weight <- edgelists.nozeros[[i]]$edgeweights
	if(length(isolated.nodes) >= 1){
		static.graph.list[[i]] <- static.graph <- add.vertices(static.graph.weighted, nv =
																 length(isolated.nodes), name =
																 as.character(isolated.nodes))
		} else {
				static.graph.list[[i]] <- static.graph <- static.graph.weighted
			}
		}
	}
}

#-- write out static graph list --#
write.path <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/GraphObjectData_27Sept2013/StaticGraphList"

write.path.edgelist <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/GraphObjectData_27Sept2013/EdgelistsNozeros"

dput(static.graph.list, paste(write.path))
dput(edgelists.nozeros, paste(write.path.edgelist))

#-- let j = 8 3-week timesteps, starting at day = 0 --#
#-- let i = number of nodes --#
#max.ts <- 8
##-- specify edgelist of interest --#
#still.alive <- edgelist.ofinterest <- el.ofinterest <- chronic.nodes  <- vector("list", length = 8)
#edgelist.ofinterest[[1]] <- edgelists.nozeros[[1]]
#el.ofinterest[[1]] <- cbind(edgelist.ofinterest[[1]]$Ind1,
#														edgelist.ofinterest[[1]]$Ind2)
#

#-- function to transmit disease through graph --#
disease.through.graph <- function(static.graph, prop.chronic,
																	P.mort.given.chronic.contact, max.ts,
																	edgelists.nozeros){

	still.alive <- edgelist.ofinterest <- static.graph.list <- el.ofinterest <- chronic.nodes  <- vector("list", length = max.ts)
	edgelist.ofinterest[[1]] <- edgelists.nozeros[[1]]
	el.ofinterest[[1]] <- cbind(edgelist.ofinterest[[1]]$Ind1,
														edgelist.ofinterest[[1]]$Ind2)

	static.graph.list[[1]] <- static.graph
	graph.edgelist(el.ofinterest[[1]], directed = F)
	prop.chronic <- prop.chronic 
	P.mort.given.chronic.contact <- P.mort.given.chronic.contact 
	surviving <- rep(NA, max.ts)

	#-- main loop -#
	for(j in 1:max.ts){
		if(j == 1){
			#-- in first timestep, select initially infected nodes --#
			chronic.nodes[[j]] <- rbinom(length(V(static.graph.list[[1]])$name), 1, prop.chronic)
			chronic.names <- V(static.graph.list[[1]])$name[chronic.nodes[[j]] == 1]
			new.I.status <- rep(NA, length(V(static.graph.list[[j]])$name))

			for(i in 1:length(levels(factor(V(static.graph.list[[j]])$name)))){
			#-- check to see if this node is already infected --#
			#-- if this node ISN'T already infected, do this: --#
			if(! (V(static.graph.list[[j]])$name[i] %in% chronic.names)){ 
			#-- extract all edges linking to node i --#
				k <- subset(edgelist.ofinterest[[j]], as.character(Ind1) ==
									V(static.graph)$name[i] | as.character(Ind2) ==
									V(static.graph)$name[i])
				#-- label all nodes with their infection status --#
				m <- subset(k, as.character(Ind1) %in% chronic.names | as.character(Ind2)
									%in% chronic.names)
				#-- sum edgeweights to infected nodes, multiply by
				#-- P.mort.given.chronic.contact --#
				prob.become.infected <- sum(m$edgeweights) * P.mort.given.chronic.contact
				infection.trial <- rbinom(1, 1, prob.become.infected)
				new.I.status[i] <- ifelse(infection.trial == 1, "I", "S")
			} else {
			new.I.status[i] <- "Dead"
			}
		}
		still.alive[[j]] <- subset(data.frame(cbind(V(static.graph)$name, new.I.status)),
																	 new.I.status != "Dead")
	names(still.alive[[j]]) <- c("Ind", "Status")
	} #-- end if(j = 1) == T) --# 
		else { #-- j > 1 == T --#
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
																 as.character(still.alive[[j - 1]]$Status) ==
																 "I")$Ind 
			for(i in 1:length(levels(factor(V(static.graph.list[[j]])$name)))){
			#-- check to see if this node is already infected --#
			#-- if this node ISN'T already infected, do this: --#
				if(! (levels(factor(V(static.graph.list[[j]])$name))[i] %in% chronic.names)){
					#-- extract all edges linking to node i --#
					k <- subset(edgelist.ofinterest[[j - 1]], as.character(Ind1) ==
									V(static.graph)$name[i] | as.character(Ind2) ==
									V(static.graph)$name[i])
				#-- label all nodes with their infection status --#
					m <- subset(k, as.character(Ind1) %in% chronic.names | as.character(Ind2)
								%in% chronic.names)
				#-- sum edgeweights to infected nodes, multiply by
				#-- P.mort.given.chronic.contact --#
					prob.become.infected <- sum(m$edgeweights) * P.mort.given.chronic.contact
					infection.trial <- rbinom(1, 1, prob.become.infected)
					new.I.status[i] <- ifelse(infection.trial == 1, "I", "S")
				} else {
						new.I.status[i] <- "Dead"
					}
			}
			still.alive[[j]] <-
			subset(data.frame(cbind(levels(factor(edgelist.ofinterest[[j]]$Ind1)), new.I.status)),
																	 new.I.status != "Dead")
			names(still.alive[[j]]) <- c("Ind", "Status")
			}	
		}
		#-- update graph by removing all Dead nodes --#
		edgelist.ofinterest[[j + 1]] <- subset(edgelist.ofinterest[[j]], as.character(Ind1)
													 %in% as.character(still.alive[[j]]$Ind) & as.character(Ind2) %in%
													 as.character(still.alive[[j]]$Ind))
		el <- cbind(as.character(edgelist.ofinterest[[j + 1]]$Ind1),
							as.character(edgelist.ofinterest[[j + 1]]$Ind2))
		static.graph.list[[j + 1]] <- graph.edgelist(el, directed = F)
		E(static.graph.list[[j + 1]])$weight <- edgelist.ofinterest[[j + 1]]$edgeweights
#	}
		surviving[j] <- length(levels(factor(V(static.graph.list[[j]])$name)))
	}
return(surviving)
}

#-- replicate disease.through.graph many times. --#
batch.fun <- function(static.graph = static.graph, prop.chronic, P.mort.given.chronic.contact, max.ts,
											edgelists.nozeros, reps){
	lamb.mort.vec <- test <- tot.ewe.vec <- rep(NA, reps)
	for(r in 1:reps){
#		lamb.mort.vec[r] <- disease.through.graph(static.graph = static.graph, 
#																							prop.chronic = prop.chronic,
#																							P.mort.given.chronic.contact =
#																							P.mort.given.chronic.contact, 
#																							max.ts = max.ts, 
#																							edgelists.nozeros =
#																							edgelists.nozeros)[max.ts]
		test[r] <- try(disease.through.graph(static.graph = static.graph, 
																							prop.chronic = prop.chronic,
																							P.mort.given.chronic.contact =
																							P.mort.given.chronic.contact, 
																							max.ts = max.ts, 
																							edgelists.nozeros =
																							edgelists.nozeros)[max.ts])
		lamb.mort.vec[r] <- ifelse(class(test[r]) != "try-error", test[r], NA)
		tot.ewe.vec[r] <- length(levels(factor(V(static.graph)$name))) 
	print(r)
	}
	out.mat <- cbind(lamb.mort.vec, tot.ewe.vec)
	return(out.mat)
}

min.Pmort <- .2  
max.Pmort <- .9
mortstep <- .1
min.Pchron <- .2
max.Pchron <- .9
chronstep <- .1
static.graph <- static.graph
max.ts <- 8
reps <- 10
edgelists.nozeros <- edgelists.nozeros 

test1 <- disease.through.graph(static.graph, prop.chronic = .2, P.mort.given.chronic.contact = .8, max.ts = 8, edgelists.nozeros) 

batch.test <- batch.fun(static.graph, prop.chronic = .2,
												P.mort.given.chronic.contact = .8, max.ts = 8,
												edgelists.nozeros, reps = 100) 

##-- brute force!! --#
#mort.1.chron.1 <- batch.fun(static.graph, prop.chronic = .1,
#												P.mort.given.chronic.contact = .1, max.ts = 8,
#												edgelists.nozeros, reps = 100) 
#
#mort.2.chron.1 <- batch.fun(static.graph, prop.chronic = .1,
#												P.mort.given.chronic.contact = .2, max.ts = 8,
#												edgelists.nozeros, reps = 100) 
#
#mort.3.chron.1 <- batch.fun(static.graph, prop.chronic = .1,
#												P.mort.given.chronic.contact = .3, max.ts = 8,
#												edgelists.nozeros, reps = 100) 
#
#mort.4.chron.1 <- batch.fun(static.graph, prop.chronic = .1,
#												P.mort.given.chronic.contact = .4, max.ts = 8,
#												edgelists.nozeros, reps = 100) 
#
#mort.4.chron.1.b2 <- batch.fun(static.graph, prop.chronic = .1,
#												P.mort.given.chronic.contact = .4, max.ts = 8,
#												edgelists.nozeros, reps = 100) 
#
#mort.5.chron.1 <- batch.fun(static.graph, prop.chronic = .1,
#												P.mort.given.chronic.contact = .5, max.ts = 8,
#												edgelists.nozeros, reps = 100) 
#
#mort.6.chron.1 <- batch.fun(static.graph, prop.chronic = .1,
#												P.mort.given.chronic.contact = .6, max.ts = 8,
#												edgelists.nozeros, reps = 100) 
#
#mort.7.chron.1 <- batch.fun(static.graph, prop.chronic = .1,
#												P.mort.given.chronic.contact = .7, max.ts = 8,
#												edgelists.nozeros, reps = 100) 
#
#mort.7.chron.1.b2 <- batch.fun(static.graph, prop.chronic = .1,
#												P.mort.given.chronic.contact = .7, max.ts = 8,
#												edgelists.nozeros, reps = 100) 
#
#mort.8.chron.1 <- batch.fun(static.graph, prop.chronic = .1,
#												P.mort.given.chronic.contact = .8, max.ts = 8,
#												edgelists.nozeros, reps = 100) 
#
#mort.9.chron.1 <- batch.fun(static.graph, prop.chronic = .1,
#												P.mort.given.chronic.contact = .9, max.ts = 8,
#												edgelists.nozeros, reps = 100) 
#
#mort.9.chron.1.b2 <- batch.fun(static.graph, prop.chronic = .1,
#												P.mort.given.chronic.contact = .9, max.ts = 8,
#												edgelists.nozeros, reps = 100) 
#
#mort.9.chron.1.b3 <- batch.fun(static.graph, prop.chronic = .1,
#												P.mort.given.chronic.contact = .9, max.ts = 8,
#												edgelists.nozeros, reps = 100) 

param.space.fun <- function(min.Pmort, max.Pmort, mortstep, min.Pchron, max.Pchron, chronstep, static.graph, max.ts, reps, edgelists.nozeros){ 
	mort.steps <- seq(min.Pmort, max.Pmort, by = mortstep) 
	chron.steps <- seq(min.Pchron, max.Pchron, by = chronstep) 
	param.mat <- expand.grid(mort.steps, chron.steps) 
	out.list <- vector("list", dim(param.mat)[1]) 
	for(p in 1:dim(param.mat)[1]){ 
		#-- need to build a storage object... --# 
		out.list[[p]] <- batch.fun(static.graph, prop.chronic = param.mat[p, 1], 
															 P.mort.given.chronic.contact = param.mat[p, 2], 
															 max.ts = max.ts, edgelists.nozeros = edgelists.nozeros, 
															 reps = reps) 
	} 
	return.list <- list(param.mat = param.mat, out.list = out.list) 
	return(return.list) 
} 

test1 <- disease.through.graph(static.graph, prop.chronic = .2, P.mort.given.chronic.contact = .8, max.ts = 8, edgelists.nozeros) 

batch.test <- batch.fun(static.graph, prop.chronic = .2,
												P.mort.given.chronic.contact = .8, max.ts = 8,
												edgelists.nozeros, reps = 100) 

param.space.test <- param.space.fun(min.Pmort = .2, max.Pmort = .9, mortstep =
																		.1, min.Pchron = .2, max.Pchron = 1,
																		chronstep = .1, static.graph =
																		static.graph, max.ts = 8, reps = 150, edgelists.nozeros = edgelists.nozeros)

#-- to process param.space.test --#
	#-- 1) loop through out.list; calculate prop successes in each run --#
	#-- 2) determine proportion of runs falling within 1 sd of observed repro rate 
		#-- (for BB in 1997, this is 5/12; SE(p-hat) = \sqrt{p(1-1)/n} = 0.1423 --#
		#-- so acceptable values fall between 5/12 - 0.1423 = .2744 & 5/12 + 0.1423
		#-- = 0.5590.  Round to more extreme full-interger values, which are 3 and
		#-- 7 here.  Any outcome between 3 and 7 is deemed "consistent" with the
		#-- empirical patterns --#

param.mat.out <- as.data.frame(param.space.test$param.mat)
names(param.mat.out) <- c("mort.rate", "chron.rate")
param.mat.out$sims.consistent <- rep(NA, dim(param.mat.out)[1])
for(i in 1:dim(param.space.test$param.mat)[1]){
	k <- table(as.numeric(as.character(param.space.test$out.list[[i]])) %in%
						 3:7)["TRUE"]
	param.mat.out$sims.consistent[i] <- ifelse(is.na(k), 0, k / 150)
}

#-- pull range of chron.rate that exceeds 60% consistency --#
chron.range <- range(subset(param.mat.out, sims.consistent >= .75)$chron.rate)
mort.range <- range(subset(param.mat.out, sims.consistent >= .75)$mort.rate)

#-- plot param.mat.out$sims.consistent as heat map --#
image.mort <- as.numeric(as.character(levels(factor(param.space.test$param.mat[
																								 ,1]))))
image.chron <- as.numeric(as.character(levels(factor(param.space.test$param.mat[
																								 ,2]))))
image.simsconsistent <- matrix(param.mat.out$sims.consistent, nrow = length(image.chron),
																	ncol = length(image.mort), byrow = T)

require(graphics)
image(x = image.mort, y = image.chron,
			z = image.simsconsistent)

require(ggplot2)
p <- ggplot(param.mat.out, aes(x = chron.rate, y = mort.rate))
plot.path <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Plots/ComponentPlots/NetworksByPopYear/RevisedPlots_18Sept2013/BlackButte/HeatMaps/BB97.jpg"
full.plot <- p + geom_tile(aes(fill = sims.consistent)) + scale_fill_gradient(low = "white",
																																 high =
																																 "black") +
theme_bw()

ggsave(paste(plot.path, sep = ""), full.plot)
