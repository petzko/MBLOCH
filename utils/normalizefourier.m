function nconst = normalizefourier(NFFT,win,transform)

win = reshape(win,length(win),1);
Y = transform(ones(length(win),1).*win,NFFT);  
nconst = abs(Y(1));

end