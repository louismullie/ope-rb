module OPE
  
  class PRNG
  
    require 'bigdecimal'
    require 'drbg-rb'
    
    def initialize(coin)
      coin = coin.to_s
      coin += '0' while coin.bytesize < 24
      @prng = DRBG::HMAC.new(coin, 128)
    end
    
    def draw
      n = @prng.generate(4, 128).unpack('H*')[0]
      n.byteslice(0, 8).hex.to_f / 2**32
    end
    
  end
  
end