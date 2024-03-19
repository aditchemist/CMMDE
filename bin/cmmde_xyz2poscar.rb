#!/usr/bin/env ruby

require 'ostruct'

input_file = ARGV[0]

def parsexyz(infile)

    res = OpenStruct.new
    res.coords = []
    res.symbols = []
    res.natoms = 0
    res.title = "POSCAR"
    res.tv = []
    res.succ = true
    res.nats = []

    counting_map = {}
    counting_map.default = 0

    natoms = infile.gets.split[0].to_i
    title = infile.gets
    natoms.times do
        line = infile.gets
        arr = line.split
        res.symbols << arr[0]
        res.coords << arr
        counting_map[arr[0]] += 1
    end 
    3.times do
        if line = infile.gets and line.include?("TV")
            res.tv << line.split(" ")[1..-1]
        else
            res.tv = [ [20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0] ]
            break
        end
    end
    res.symbols.uniq!
    res.symbols.each do |sym|
        res.nats << counting_map[sym]
    end
    return res
end


def writePOSCAR(res)
    if (res.title.size == 0)
        puts "POSCAR"
    else
        puts res.title
    end
    puts 1.0
    (0..2).each do |i|
        puts res.tv[i].map{|x| "%16.12f" % x}.join(" ")
    end
    puts res.symbols.join(" ")
    puts res.nats.join(" ")
    puts "Cart"
    res.coords.each do |line|
        puts "#{line[1..-1].map{|x| "%16.12f" % x}.join(" ")}   #{line[0]}"
    end
end




File.open(input_file) do |infile|
    num = 0
    last_obj = OpenStruct.new
    while infile.eof? == false
        obj = parsexyz(infile)
        num += 1
        if obj.succ
            last_obj = obj
#            last_obj.marshal_load(obj.marshal_dump())
        end
    end
    writePOSCAR(last_obj)
end
exit


File.open(output_file, "w") do |outfile|
    
    outfile.puts [natoms, "C"].join(" ")
    outfile.puts symbols.join(" ")
    lines.each_with_index do |line, index|
        outfile.puts [index+1, symbols.index(line[0])+1, line[1..3]].join("  ")
    end
end
