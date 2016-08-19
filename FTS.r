path <- "/home/petronio/dados/Dropbox/Doutorado/Disciplinas/AdvancedFuzzyTimeSeriesModels/rfts/"
files <- c("common.r","partitioner.r","ftscommon.r","song.r",
"chen.r","yu.r","ismail_efendi.r","sadaei.r",
"benchmarks.r","ifts.r","hofts.r","sfts.r")

for (i in files) {
	source(paste(path,i,sep=""))
}
        

