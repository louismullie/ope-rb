module OPE
  
  class PRNG
  
    require 'bigdecimal'
    require 'openssl'
    
    def initialize(coin)
      @prng = Random.new(coin)
    end
    
    # Returns a random number in [0, 1]
    def draw
      @prng.rand
    end
    
  end
  
end