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

=begin
File.open(ARGV[1], "r") do |ref_file|
	reference = Bio::FlatFile.auto(ref_file)
	reference.each{|each_chr|
		chr.push([each_chr.entry_id, each_chr.seq.to_str])
	}
end
=end

bp_list = Hash.new {|h,k| h[k] = []}

raw_sv_list.each{|each_sv|
	bp_list.store(each_sv[0], bp_list[each_sv[0]].push(each_sv[1].to_i))
	bp_list.store(each_sv[3], bp_list[each_sv[3]].push(each_sv[4].to_i))
=begin
	if each_sv[7] == "DEL"
		seq = "#{each_sv[0]}:#{each_sv[1].to_i}-#{each_sv[4].to_i}"
		segment = `samtools faidx #{ARGV[1]} #{seq}`
		#segment.split("\n").drop(1).join("").upcase.length
	end
	if each_sv[7] == "INS"
		seq = "#{each_sv[0]}:#{each_sv[1].to_i}-#{each_sv[4].to_i}"
		segment = `samtools faidx #{ARGV[1]} #{seq}`
		#segment.split("\n").drop(1).join("").upcase.length
	end
=end
}
#bp_list contains chr and locus
#{chr, [###,###,###,###,...,###]}
uniq_id = 0
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
		#puts segment.split("\n").drop(1).join("").upcase
		puts "S\t#{uniq_id}\t#{segment.split("\n").drop(1).join("").upcase}"
		uniq_id = uniq_id + 1
		puts "L\t#{i}\t+\t#{i + 1}\t+\t0M"
		
		#p segment.split("\n").drop(1).join("").upcase
		i = i + 1
	end
	#p list[1][list[1].length]
}
bp_list.each{|list|

}
