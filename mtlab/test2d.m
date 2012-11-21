 close all
[U,V,Y,W,R1,R2,Is]=fol2d3;

  m1=makemovie2d(20,20,U,V,Is,'bmp','wnt');
  m2=makemovie2d(20,20,W,Y,Is,'Noggin','Dkk');
  movie2avi(m1,'ligands.avi')
  movie2avi(m2,'antagonists.avi')