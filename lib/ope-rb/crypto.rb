module OPE
  
  module Crypto

    def hash_concat(*args)
      
      lens = args.map { |arg| [1].pack("L") }
      hash((lens + args).join(''))
      
    end
    
    def hash(data)
      
      OpenSSL::Digest::SHA256.digest(data)
      
    end
    
    def create_cipher
      
      @cipher = OpenSSL::Cipher::AES.new(128, :ECB)
      @cipher.encrypt; @cipher.key = @key
      
    end
    
    def cipher_encrypt(data)
      
      @cipher.update(data) + @cipher.final
      
    end
  
  end
  
end