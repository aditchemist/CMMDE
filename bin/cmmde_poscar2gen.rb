#!/usr/bin/env ruby


require 'matrix'

lines = []
while line = gets
    lines << line
end  


title = lines.shift

factor = lines.shift.to_f
basis_vectors = []
3.times do
   basis_vectors << lines.shift
end

basis = basis_vectors.map{|x| x.split.map{|y| y.to_f}}

typenames = lines.shift.split
nats = lines.shift.split.map{|x|x.to_i}


selective_dynamics = false

if (lines[0].downcase().start_with?("selective"))
    selective_dynamics = true
    lines.shift
end

mode_str = lines.shift
mode_type = "S"
if mode_str[0..0].upcase == "D"
    mode_name = "Relative"
    mode_type = "F"
else
    mode_name = "Angstrom"
end

coords = []
nats.each_with_index do |n, i|
    n.times do 
        line = lines.shift.split(" ")
        coords << [i+1, line[0..2]]
    end
end


final = basis.map{|line| (Vector.elements(line)*factor).to_a}

nat = nats.inject(0){|sum,x| sum + x }

puts "#{nat} #{mode_type}"
puts "#{typenames.join(" ")}"
coords.each_with_index do |line, j|
    puts "#{j+1} #{line.join("  ")}"
end

puts "0.0 0.0 0.0"
final.each do |line|
    puts line.map{|x| "%20.14f" % x}.join("  ")
end
exit

puts "
    TypeNames = { #{typenames.join(" ") } }
    TypesAndCoordinates [#{mode}] = {
"
    coords.each do |line|
        puts "        #{line.join("  ")}"
    end
puts "    }
    Periodic = Yes
    LatticeVectors [Angstrom] = {
"
puts "    }
"


