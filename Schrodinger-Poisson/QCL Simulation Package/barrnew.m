 layers = [48  76 23 74 40 162 33 92 ];
 n = length(layers)
 
 layersnew =layers;
 
 hold = 0.1126; 
 hnew = 0.52; 
 
 ratio = hold/hnew; 
 
 for i =1:n
    layersnew(i) = ratio*layers(i); 
 end
 
 barriers = layers(1:2:n);
 layersnew = round(layersnew); 
 barriersnew = layersnew(1:2:n);