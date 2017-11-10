#!/usr/bin/env ruby
require 'bio'
require 'csv'
if ARGV.size < 2 then
	puts "pcf2ggf.rb <xxx.pcf> <refarence.fa>"
	exit(1)
end
bp = []
CSV.foreach(ARGV[0]) do |row|
	bp.push(row)
end
p bp

