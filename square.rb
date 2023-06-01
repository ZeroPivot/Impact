# wolf moon bolt, wolf+moon+bolt

def num_num(char)
  letters = {"a" => 1,
             "b" => 2,
             "c" => 3,
             "d" => 4,
             "e" => 5,
             "f" => 6,
             "g" => 7,
             "h" => 8,
             "i" => 9,
             "j" => 1,
             "k" => 2,
             "l" => 3,
             "m" => 4,
             "n" => 5,
             "o" => 6,
             "p" => 7,
             "q" => 8,
             "r" => 9,
             "s" => 1,
             "t" => 2,
             "u" => 3,
             "v" => 4,
             "w" => 5,
             "x" => 6,
             "y" => 7,
             "z" => 8,}
  letters[char]
  end


def mczv(🐺, 🌕, ⚡, 🐺🌕⚡)
  
  a = 🐺
  b = 🌕
  c = ⚡
  d = a+b+c if !🐺🌕⚡
  d = 🐺🌕⚡
  e = c-2
  f = d+2
  g = a-2
  h = b+2
  i = d+1
  j = c+1
  k = b-1
  l = a-1
  m = b+1
  n = a-3
  o = d+3
  p = c-1

  puts "#{a} #{b} #{c} #{d}"
  puts "#{e} #{f} #{g} #{h}"
  puts "#{i} #{j} #{k} #{l}"
  puts "#{m} #{n} #{o} #{p}"

  
end

loop do

  puts "(1) Letters or (2) numbers?"
  answer = gets.chomp
  if answer == "1"
    puts "Enter the 4 letters"
    split = gets.chomp.split("")
    a = num_num(split[0])  
    b = num_num(split[1])    
    c = num_num(split[2])  
    d = num_num(split[3])
    puts "#{split.join("")}"
    mczv(a, b, c, d)
  end
  if answer == "2"
    puts "Enter the first number"
    a = gets.chomp.to_i
    puts "Enter the second number"
    b = gets.chomp.to_i
    puts "Enter the third number"
    c = gets.chomp.to_i
    puts "Enter the fourth number"
    d = gets.chomp.to_i
    mczv(a, b, c, d)
  end
end


# Q: how so I configure user and email in git?
# A: git config --global user.name "Your Name"