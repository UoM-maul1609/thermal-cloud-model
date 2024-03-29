function pi=svp(T,flag,flag2)
%usage: SVPIce(T,flag,flag2);
% flag='marti'
% flag='teten'
% flag2='ice'
% flag2='liq'


switch lower(flag2)
case 'ice'
    switch lower(flag)
    case 'buck2'
        pi = 100.*6.1115 .* exp((23.036 - (T-273.15)./ 333.7) .*(T-273.15) ./ (279.82 + (T-273.15)))  ;
    case 'buck'
        pi = 100.*6.1115 .* exp(22.452 .* (T-273.15) ./ (272.55+(T-273)));  

    case 'goff'
        pi =  100.*10.^(-9.09718.* (273.16./T - 1)       ...                                 
                   - 3.56654 .*log10(273.16./ T) ...
                   + 0.876793 .*(1 - T./ 273.16) ...
                   + log10(6.1071) );
        
    case 'marti'
        pi = 10.^((-2663.5 ./ T) + 12.537 );   
    case 'teten'
        pi = 100.*10.^(9.5 .*(T-273.15) ./ (T-273.15+265.5) + 0.7858  ) ;  
    case 'hyland'
        pi =  exp(-0.56745359e4 ./ T     ...                                            
              + 0.63925247e1 ...
              - 0.96778430e-2 .*T ...
             + 0.62215701e-6 .*T.^2 ...
             + 0.20747825e-8 .*T.^3 ...
             - 0.94840240e-12 .*T.^4 ...
             + 0.41635019e1 .*log(T) );
    case 'murphy'
        pi = exp(9.554605 - 5722.796./T + 3.5291623.*log(T) - 0.00727374.*T);
    case 'magnus'
        pi=610.7.*exp((22.44.*T-6.1186e3)./(T-0.75));
    case 'clausius'
        pi=611.73.*exp(2.501e6./461.5.*(1./273.16-1./T));
    otherwise
        error('invalid');
    end
case 'liq'
    switch lower(flag)
    case 'goff'
        pi =  100.*10.^(-7.90298 .*(373.16./T-1)    ...                       
                    + 5.02808 .*log10(373.16./T) ...
                    - 1.3816e-7 .*(10.^(11.344 .*(1-T./373.16))  -1) ...
                   + 8.1328e-3 .*(10.^(-3.49149 .*(373.16./T-1))  -1) ...
                   + log10(1013.246) );
        
    case 'bolton'
        pi = 100.*6.112 .*exp(17.67 .* (T-273.15) ./ (T-273.15+243.5));
    case 'roger'
        pi = 2.53e11 * exp(-5.42e3./(T));
    case 'buck2'
        pi = 100.*6.1121  .*exp((18.678 - (T-273.15)./ 234.5).* (T-273.15) ./ (257.14 + (T-273.15)));
    case 'buck1'
        pi = 100.*6.1121 .*exp(17.502 .*(T-273.15)./ (240.97 + T-273.15));

    case 'wmo'
        pi = 100.*10.^( 10.79574 .*(1-273.16./T)        ...                       
                    - 5.02800 .*log10(T./273.16) ...
                    + 1.50475e-4 .*(1 - 10.*(-8.2969.*(T./273.16-1))) ...
                    + 0.42873e-3 .*(10.*(+4.76955.*(1-273.16./T)) - 1) ...
                    + 0.78614 );
    case 'hyland'
        pi =  exp(-0.58002206e4 ./ T  ...                                     
              + 0.13914993e1 ...
              - 0.48640239e-1 .* T ...
              + 0.41764768e-4 .* T.^2 ...
              - 0.14452093e-7 .* T.^3 ...
              + 0.65459673e1 .* log(T)); 
            
    case 'sonntag'
        pi =  100.*exp(-6096.9385 ./ T  ...                        
                 + 16.635794 ...
                 - 2.711193e-2 .* T ...
                 + 1.673952e-5 .* T.^2 ... 
                 + 2.433502 .* log(T)); 
    case 'teten'
        pi = 100.*10.^(7.5 .*(T-273.15) ./ (T-273.15+237.3) + 0.7858  ) ;  
    case 'clausius'
        pi=611.73.*exp(2.834e6./461.5.*(1./273.16-1./T));
    case 'magnus'
        pi=610.7.*exp((17.38.*T-4.7473e3)./(T-34.15));
    otherwise
        error('invalid');
    end
end    
