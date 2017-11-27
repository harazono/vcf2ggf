#!/usr/bin/env ruby
require 'bio'
require 'csv'

if ARGV.size < 2 then
	puts "pcf2ggf.rb <xxx.pcf> <reference.fa>"
	exit(1)
end

raw_sv_list = []
CSV.foreach(ARGV[0]) do |row|
	raw_sv_list.push(row)
end


#get all break points and store in "bp_list"
#bp_list contains chr and locus
#{chr, [###,###,###,###,...,###]}
bp_list = Hash.new {|h,k| h[k] = []}
raw_sv_list.each{|each_sv|
	bp_list.store(each_sv[0], bp_list[each_sv[0]].push(each_sv[1].to_i))
	bp_list.store(each_sv[3], bp_list[each_sv[3]].push(each_sv[4].to_i))
}


#show SV's breakpoints as uniq_id(in each chromosome)
sv_list = []
raw_sv_list.each{|each_sv|
	chr_bp_list1 = bp_list[each_sv[0]].sort
	chr_bp_list2 = bp_list[each_sv[3]].sort
	sv_list.push([each_sv[0], chr_bp_list1.find_index{|n| n == each_sv[1].to_i}, each_sv[2], each_sv[3], chr_bp_list2.find_index{|n| n == each_sv[4].to_i}, each_sv[5], each_sv[6].to_i, each_sv[7], each_sv[8]])
}

# cut reference sequence, make Segments and connect Segments with Link.
# build reference graph genome
uniq_id = 0
ref_segment_name_list = ""
ref_cigar_list = ""
bp_list.each{|list|
	list[1] = list[1].sort
	pre_bp = 0
	nxt_bp = 0
	i = 0
	while list[1][i] != nil do
		pre_bp = list[1][i]
		if list[1][i + 1] != nil then
			nxt_bp = list[1][i + 1]
		else
			nxt_bp = 270000000
		end
		seq = "#{list[0]}:#{pre_bp}-#{nxt_bp}"
		segment = `samtools faidx #{ARGV[1]} #{seq}`
		puts "S\t#{uniq_id}\t#{segment.split("\n").drop(1).join("").upcase}"
		puts "L\t#{uniq_id}\t+\t#{uniq_id + 1}\t+\t0M"
		ref_segment_name_list << uniq_id.to_s << "+,"
		ref_cigar_list << "0M,"
		uniq_id = uniq_id + 1
		i = i + 1
	end
}
puts "P\tref\t#{ref_segment_name_list.chop}\t#{ref_cigar_list.chop}"

