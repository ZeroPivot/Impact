# frozen_string_literal: true
# require a directory of files that end in .rb

def require_dir(dir)
  Dir.glob(File.join(dir, '**', '*.rb')).each do |file|
    require_relative file
    puts "require_relative #{file}"
  end
end

puts "blah's require_dir.rb"


