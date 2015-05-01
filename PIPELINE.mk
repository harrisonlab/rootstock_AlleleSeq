
BASE=xxx

PL:=.
SNPS:=$(BASE)/x.snv
CNVS:=$(BASE)/x.snv.cnv 
BNDS:=hits.bed
MAPS:=$(BASE)/chr%s_x.map
FDR_SIMS:=5
FDR_CUTOFF:=0.1


sourcefiles := $(RNA_1)
countfiles := $(subst .fq.trim,.cnt,$(sourcefiles))

%.cnt:%.fq.trim 	
	bash -c "python $(PL)/MergeBowtie.py \
           <($(PL)/filter_input.sh $(PL) $(RNA_1) $(RNA_2) a1; bowtie2 --no-hd --no-unal -f -x $(BASE)/pat2/m116_paternal_index -1 a1_1 -2 a1_2|grep -v XS |cut -f 1,2,3,4,10) \
           <($(PL)/filter_input.sh $(PL) $(RNA_1) $(RNA_2) a2; bowtie2 --no-hd --no-unal -f -x $(BASE)/mat2/m116_maternal_index -1 a2_1 -2 a2_2|grep -v XS |cut -f 1,2,3,4,10) \
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
