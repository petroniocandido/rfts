path <- "/home/petronio/dados/Dropbox/Doutorado/Disciplinas/AdvancedFuzzyTimeSeriesModels/rfts/"
files <- c("common.r","partitioner.r","ftscommon.r","song.r","chen.r","benchmarks.r","pfts.r")

for (i in files) {
	source(paste(path,i,sep=""))
}
        

