#!/usr/bin/env ruby

require 'matrix'


genfile = ARGV[0]

type = "C"

symbols = []
index = []
coords = []
nat = 0

lat_vector = []

File.open(genfile, "r") do |file|
    line = file.gets
    arr = line.split
    if arr.size > 1
        type = arr[1]
    end
    nat = arr[0].to_i
    symbols = file.gets.split
    (1..nat).each do |i|
        line = file.gets
        arr = line.split
        index << arr[1].to_i
        coords << Vector.elements(arr[2..4].map{|x| x.to_f})
    end
    if ( type.upcase != "C")
        file.gets
        (1..3).each do |i|
            arr = file.gets.split.map{|x| x.to_f}
            lat_vector << Vector.elements(arr)
        end
    end
end
 

puts nat
puts 
if type.upcase == "C"
    (1..nat).each do |i|
        puts [symbols[index[i-1]-1], coords[i-1].to_a].to_a.join("  ")
    end
elsif type.upcase == "F"
    (1..nat).each do |i|
        puts [symbols[index[i-1]-1], (coords[i-1][0]*lat_vector[0]+coords[i-1][1]*lat_vector[1]+coords[i-1][2]*lat_vector[2]).to_a.map{|x| "%16.12f"%x}].join("  ")
    end
elsif type.upcase == "S"
    (1..nat).each do |i|
        puts [symbols[index[i-1]-1], coords[i-1].to_a.map{|x| "%16.12f"%x}.join(" ")].to_a.join("  ")
    end
end
if ( type.upcase != "C")
puts "TV #{lat_vector[0].to_a.map{|x| "%16.12f"%x}.join(" ")}"
puts "TV #{lat_vector[1].to_a.map{|x| "%16.12f"%x}.join(" ")}"
puts "TV #{lat_vector[2].to_a.map{|x| "%16.12f"%x}.join(" ")}"
end
