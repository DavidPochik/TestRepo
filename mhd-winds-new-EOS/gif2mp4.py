import sys
import os
import moviepy.editor as mp
gifname      = str(sys.argv[1])
gifname_full = gifname+'.gif'
clip = mp.VideoFileClip(gifname_full)
clip.write_videofile(gifname+'.mp4')
