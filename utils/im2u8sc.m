function im = im2u8sc(im, caxis)
im = uint8(256*(im/caxis(2)-caxis(1)));
