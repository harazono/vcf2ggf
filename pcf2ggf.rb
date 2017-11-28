#!/usr/bin/env ruby
#require 'bio'
require 'csv'

if ARGV.size < 2 then
	puts "pcf2ggf.rb <xxx.pcf> <reference.fa>"
	exit(1)
end

puts "H\tVN:Z:1.0"

raw_sv_list = []
CSV.foreach(ARGV[0], {headers: true}) do |row|
	raw_sv_list.push(row) if row[0].start_with?("#") == false
end

bp_list = Hash.new {|h,k| h[k] = []}
raw_sv_list.each{|each_sv|
	bp_list.store(each_sv[0], bp_list[each_sv[0]].push(each_sv[1].to_i))
	bp_list.store(each_sv[3], bp_list[each_sv[3]].push(each_sv[4].to_i))
}
bp_list.each{|list|
	list[1].sort!
}


#chr 1 <S1>bp1<S2>bp2<S3>
#chr 2 <S4>bp1<S5>bp2<S6>bp3<S7>
#get all break points and store in bp_uniq_id_list
#chrN=>locus=>uniq_id
#{"chrN"=>{locus=>uniq_id}}
bp_uniq_id_list = Hash::new()
raw_sv_list.each{|each_sv|
	key = each_sv[0]
	bp_uniq_id_list[key] = Hash::new() if bp_uniq_id_list[key] == nil
	bp_uniq_id_list[key][each_sv[1].to_i] = bp_list[each_sv[0]].find_index{|n| n == each_sv[1].to_i}
	key = each_sv[3]
	bp_uniq_id_list[key] = Hash::new() if bp_uniq_id_list[key] == nil
	bp_uniq_id_list[key][each_sv[4].to_i] = bp_list[each_sv[3]].find_index{|n| n == each_sv[4].to_i}
}


sv_list = []
raw_sv_list.each{|each_sv|
	sv_list.push([each_sv[0], bp_uniq_id_list[each_sv[0]][each_sv[1].to_i] , each_sv[2], each_sv[3], bp_uniq_id_list[each_sv[3]][each_sv[4].to_i], each_sv[5], each_sv[6].to_i, each_sv[7], each_sv[8]])
}
#[["chr20", 15568, "-", "chr20", 15569, "-", 1, "INV", nil], ["chr20", 6748, "-", "chr20", 6881, "+", 1, "DUP", nil], ["chr20", 6749, "+", "chr20", 6882, "-", 1, "DEL", nil]]






# cut reference sequence, make Segments and connect Segments with Link.
# build reference graph genome
uniq_id = 1
ref_segment_name_list = ""
ref_cigar_list = ""
bp_list.each{|list|
	pre_bp = 0
	nxt_bp = 0
	i = 0
	while nxt_bp != 270000000
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
		ref_cigar_list << "#{segment.length}M,"
		uniq_id = uniq_id + 1
		i = i + 1
	end
	puts "P\t#{list[0]}\t#{ref_segment_name_list.chop}"
	#now reference genome's Path was added!
}


#add new Link line connected with Structural Variations
del_uid = 1
ins_uid = 1
inv_uid = 1
sv_list.each{|list|
	case list[7]
	when "DEL"
		from_seg = list[1] + 1#bp_uniq_id_list[list[0]].key(list[1])
		to_seg = list[4] + 1#bp_uniq_id_list[list[3]].key(list[4])
		puts "L\t#{from_seg}\t+\t#{to_seg}\t+\t0M"
		puts "P\tDEL_#{del_uid}\t#{from_seg}+,#{to_seg}+"#need to add segment length
		del_uid = del_uid + 1
	when "INS" #need to check which side INS sequence will be inserted.right side of breakpoint or left side of break point
		puts "S\t#{uniq_id}\t#{list[8].strip}"
		uniq_id = uniq_id + 1
		from_seg = list[1] + 1#bp_uniq_id_list[list[0]].key(list[1])
		to_seg = list[4] + 1#bp_uniq_id_list[list[3]].key(list[4]) + 1
		ins_seg = uniq_id
		puts "L\t#{from_seg}\t+\t#{ins_seg}\t+\t0M"
		puts "L\t#{ins_seg}\t+\t#{to_seg}\t+\t0M"
		puts "P\tINS_#{ins_uid}\t#{from_seg}+,#{ins_seg}+,#{to_seg}+"
		ins_uid = ins_uid + 1
	when "INV"
		from_seg = list[1] + 1#bp_uniq_id_list[list[0]].key(list[1])
		to_seg = list[4] + 1#bp_uniq_id_list[list[3]].key(list[4]) + 1
		seg_name = ""
		cigar_name = ""
		tmpcnt = to_seg
		puts "L\t#{from_seg}\t+\t#{to_seg}\t-\t0M"
		seg_name << from_seg.to_s << "+,"
		cigar_name << "0M,"
		while tmpcnt - 1 > from_seg
			puts "L\t#{tmpcnt}\t-\t#{tmpcnt - 1}\t-\t0M"
			seg_name << tmpcnt.to_s << "-,"
			cigar_name << "0M,"
			tmpcnt = tmpcnt - 1
		end
		seg_name << (from_seg + 1).to_s << "-," << (to_seg + 1).to_s << "+"
		cigar_name << "0M"
		puts "L\t#{from_seg + 1}\t-\t#{to_seg + 1}\t+\t0M"
		puts "P\tINV_#{inv_uid}\t#{seg_name}\t"
		inv_uid = inv_uid + 1
	when "DUP"
		#puts "DUP"
	else
		#puts "other"
	end
}








