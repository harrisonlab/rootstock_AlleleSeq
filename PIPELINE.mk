
BASE=/home/groups/harrisonlab/project_files/rootstock_genetics/m116/allele

PL:=.
SNPS:=$(BASE)/m116.snv
CNVS:=$(BASE)/m116.snv.cnv 
BNDS:=hits.bed
MAPS:=$(BASE)/chr%s_m116.map
FDR_SIMS:=5
FDR_CUTOFF:=0.1



sourcefiles := $(wildcard *.fastq.gz)
countfiles := $(subst .fastq.gz,.cnt,$(sourcefiles))

%.cnt:%.fastq.gz
	bash -c "python $(PL)/MergeBowtie.py \
           <($(PL)/filter_input.sh $(PL) $< | bowtie --best --strata -v 2 -m 1 -f $(BASE)/pat/m116_paternal_index - ) \
           <($(PL)/filter_input.sh $(PL) $< | bowtie --best --strata -v 2 -m 1 -f $(BASE)/mat/m116_maternal_index - ) \
           $(MAPS) | python $(PL)/SnpCounts.py $(SNPS) - $(MAPS) $@"

all: interestingHets.txt

check:
	@echo "checking"
	@echo $(sourcefiles)


counts.txt: $(countfiles)
	python $(PL)/CombineSnpCounts.py 5 $(SNPS) $(BNDS) $(CNVS) counts.txt counts.log $(countfiles)

# calculate false discovery rates
FDR.txt: counts.txt
	python $(PL)/FalsePos.py counts.txt $(FDR_SIMS) $(FDR_CUTOFF) > FDR.txt

interestingHets.txt: counts.txt FDR.txt
	awk -f $(PL)/finalFilter.awk thresh=$(shell awk 'END {print $$6}' FDR.txt) < counts.txt > interestingHets.txt

clean:
#	@rm -f FDR.txt interestingHets.txt counts.txt

cleanall: clean
#	@rm -f *.cnt

#.DELETE_ON_ERROR:
