module OPE
  
  module Errors
    
    class KeySizeError < StandardError
    
      def initialize
        super('Incorrect key size')
      end
    
    end
  
    class HGDParameterError < StandardError
    
      def initialize
        super('Corrupt HGD parameters')
      end
    
    end
  
    class NumBitsError < StandardError
    
      def initialize
        super('Desired no bits not a multiple of 8')
      end
    
    end
  
    class DesiredBytesError < StandardError
    
      def initialize
        super('Only 64-bit output is supported for now')
      end
    
    end
  
    class IncorrectMValueError < StandardError
    
      def initialize
        super('M should not be less than or equal to 0')
      end
    
    end
  
  end
  
end