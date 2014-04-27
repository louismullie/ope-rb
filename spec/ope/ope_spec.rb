require 'spec_helper'
require 'ope'

describe OPE do
  
  it "should work the HGD" do
    
    # puts "test: " + OPE::HGD.afc(4343423432).to_s
    # puts "test2: " + OPE::HGD.afc_native(4343423432).to_s
    
    Benchmark.bm { |x| x.report {
      
      RubyProf.start
      
      a = nil
      
      1000.times do 
      
        a = OPE::HGD.rhyper(5070602400912917605986812821503, 549756130608.0, 10141204801825835000000000000000.0, 15109037923498441947, 10)
    
      end
      
      puts a.inspect
      
      result = RubyProf.stop

      # Print a flat profile to text
      printer = RubyProf::FlatPrinter.new(result)
      printer.print(STDOUT)
      
    
     }}
     
    # OPE::HGD.rhyper_native('5070602400912917605986812821503', '549756130608', '10141204801825835000000000000000', '15109037923498441947', 10)
    
  end
=begin
  Key = ([0] * 16).pack('C*')
  
  it "should encrypt/decrypt roundtrips" do
    
    cipher = OPE::Cipher.new(Key)

    100.times do
      
      num = (rand * 20_000).floor
      
      enc_num = cipher.encrypt(num)
      dec_num = cipher.decrypt(enc_num)
      
      dec_num.should eql num
      
    end
    
  end
  
  it "should produce deterministic output" do
    
    cipher = OPE::Cipher.new(Key)

    100.times do
      
      num = (rand * 20_000).floor
      
      a = cipher.encrypt(num)
      b = cipher.encrypt(num)
      
      a.should eql b
      
    end

  end

  it "should maintain proper order" do
    
    cipher = OPE::Cipher.new(Key)
    
    plaintexts = []
    encrypted = []
    
    100.times do
      num = (rand * 100_000).floor
      plaintexts << num
      encrypted << cipher.encrypt(num)
    end
    
    while !plaintexts.empty?
      
      max_dec_pos = plaintexts.index(plaintexts.max)
      max_enc_pos = plaintexts.index(plaintexts.max)
      
      max_dec_pos.should eql max_enc_pos
      
      plaintexts.delete_at(max_dec_pos)
      encrypted.delete_at(max_enc_pos)
      
    end
      
  end
=end
end