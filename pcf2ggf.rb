#!/usr/bin/env ruby
require 'csv'

if ARGV.size < 2 then
	puts "pcf2ggf.rb <xxx.pcf> <reference.fa>"
	exit(1)
end

puts "H\tVN:Z:1.0"

#Skip header line start with "#"
raw_sv_list = []
CSV.foreach(ARGV[0], {headers: true}) do |row|
	raw_sv_list.push(row) if row[0].start_with?("#") == false
end
#raw_sv_list contains each rows.
#155,chrX,144341561,+,chrX,144341593,-,10,DEL,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#raw_sv_list sample
#[#<CSV::Row "uniq_id":"0" "source_id":"chrX" "source_breakpoint":"141906130" "source_strand":"+" "target_id":"chrX" "target_breakpoint":"141906130" "target_strand":"-" "priority":"5" "svtype":"INS" "sequence":"GAA...G">, #<CSV::Row "uniq_id":"1" "source_id":"chr20" "source_breakpoint":"35226715" "source_strand":"+" "target_id":"chr20" "target_breakpoint":"35226715" "target_strand":"-" "priority":"3" "svtype":"INS" "sequence":"ATG...TCA">, #<CSV::Row "uniq_id":"2" "source_id":"chr1" "source_breakpoint":"12762482" "source_strand":"+" "target_id":"chr1" "target_breakpoint":"12762482" "target_strand":"-" "priority":"3" "svtype":"INS" "sequence":"AAA...TTC">]

bp_list = Hash.new {|h,k| h[k] = []}
raw_sv_list.each{|each_sv|
	bp_list.store(each_sv[1], bp_list[each_sv[1]].push(each_sv[2].to_i))
	bp_list.store(each_sv[4], bp_list[each_sv[4]].push(each_sv[5].to_i)) if each_sv[8] != "INS"
}
bp_list.each{|list|
	list[1].sort!
}
#bp_list sample
#{"chrX"=>[141906130, 141906130], "chr20"=>[35226715, 35226715], "chr1"=>[12762482, 12762482]}


#chr 1 <S1>bp1<S2>bp2<S3>
#chr 2 <S4>bp1<S5>bp2<S6>bp3<S7>
#get all break points and store in bp_uniq_id_list
#chrN=>locus=>uniq_id
#{"chrN"=>{locus=>uniq_id}}
bp_uniq_id_list = Hash::new()
raw_sv_list.each{|each_sv|
	key = each_sv[1]
	bp_uniq_id_list[key] = Hash::new() if bp_uniq_id_list[key] == nil
	bp_uniq_id_list[key][each_sv[2].to_i] = bp_list[each_sv[1]].find_index{|n| n == each_sv[2].to_i}
	key = each_sv[4]
	bp_uniq_id_list[key] = Hash::new() if bp_uniq_id_list[key] == nil
	bp_uniq_id_list[key][each_sv[5].to_i] = bp_list[each_sv[4]].find_index{|n| n == each_sv[5].to_i}
}
#bp_uniq_id_list sample
#{"chr11"=>{108715009=>0}, "chr13"=>{21168402=>0}, "chr12"=>{99376451=>1, 99376310=>0}, "chr9"=>{71741204=>1, 70623503=>0}, "chr21"=>{10415117=>0}, "chr7"=>{152401907=>0}}


sv_list = []
raw_sv_list.each{|each_sv|
	sv_list.push([each_sv[0].to_i, each_sv[1], bp_uniq_id_list[each_sv[1]][each_sv[2].to_i] , each_sv[3], each_sv[4], bp_uniq_id_list[each_sv[4]][each_sv[5].to_i], each_sv[6], each_sv[7].to_i, each_sv[8], each_sv[9]])
}

# sv_list sample
#[[0, "chr11", 0, "-", "chr13", 0, "+", 6, "TRA", nil], [1, "chr12", 1, "-", "chr9", 1, "-", 4, "TRA", nil], [2, "chr12", 0, "+", "chr9", 0, "-", 3, "TRA", nil], [3, "chr21", 0, "-", "chr7", 0, "-", 4, "TRA", nil]]

# cut reference sequence, make Segments and connect Segments with Link.
# build reference graph genome

MAX = 270000000
uniq_id = 0
ref_segment_name_list = ""
segment = ""
ref_cigar_list = ""
segment_metadata = []
bp_list.each{|list|
	pre_bp = 0
	nxt_bp = 0
	i = 0
	while nxt_bp != MAX
		if i != 0 then
			pre_bp = nxt_bp
			if list[1][i] != nil then
				nxt_bp = list[1][i]
			else
				nxt_bp = MAX
			end
		else
			pre_bp = 0
			nxt_bp = list[1][i]
		end
		seq = "#{list[0]}:#{pre_bp}-#{nxt_bp}"
		segment = `samtools faidx #{ARGV[1]} #{seq}`
		puts "S\t#{uniq_id}\t"#{segment.split("\n").drop(1).join("").upcase}"
		puts "L\t#{uniq_id}\t+\t#{uniq_id + 1}\t+\t0M"
		segment_metadata.push([uniq_id, list[0], i, pre_bp, nxt_bp])

		ref_segment_name_list << uniq_id.to_s << "+,"
		ref_cigar_list << "#{segment.length}M,"
		uniq_id = uniq_id + 1
		i = i + 1
	end
	puts "P\t#{list[0]}\t#{ref_segment_name_list.chop}"
	ref_segment_name_list = ""
	#now reference genome's Path was added!
}
segment_metadata.each{|pigy|
	p pigy
}
sv_list.each{|zura|
	p zura
}
#add new Link line connected with Structural Variations
del_uid = 1
ins_uid = 1
inv_uid = 1
sv_list.each{|list|
	case list[8]#svtype
	when "DEL"
		from_seg = segment_metadata.select{|item| item[1]==list[1] && item[2] == list[2]}[0][0]
		to_seg = segment_metadata.select{|item| item[1]==list[4] && item[2] == list[5] + 1}[0][0]
		puts "L\t#{from_seg}\t+\t#{to_seg}\t+\t0M"
		puts "P\tDEL_#{del_uid}\t#{from_seg}+,#{to_seg}+"#need to add segment length
		del_uid = del_uid + 1
	when "INS" #if SV starting breakpoint is differ from SV ending point, use first one.
		if list[9] != nil
			puts "S\t#{uniq_id}\t#{list[9].strip}"#brand new segment
			uniq_id = uniq_id + 1#brand new segment's unique id
			from_seg = segment_metadata.select{|item| item[1]==list[1] && item[2] == list[2]}[0][0]
			to_seg = from_seg + 1
			ins_seg = uniq_id
			puts "L\t#{from_seg}\t+\t#{ins_seg}\t+\t0M"
			puts "L\t#{ins_seg}\t+\t#{to_seg}\t+\t0M"
			puts "P\tINS_#{ins_uid}\t#{from_seg}+,#{ins_seg}+,#{to_seg}+"
			ins_uid = ins_uid + 1
		end
	when "INV"
		from_seg = list[2] + 1#bp_uniq_id_list[list[0]].key(list[1])
		to_seg = list[5] + 2#bp_uniq_id_list[list[3]].key(list[4]) + 1
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
