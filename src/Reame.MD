Note, the frustratometer needs the pdb structure to have all the backbone and the C-beta present in order to run. If you have some of those atoms missing you will need to complete the structure. In the server we apply a simple automodel with modeller to do that, but you can use whatever other tool to do it. This server also does the work but the code is not available (at least the last time I checked) to run it locally: http://lorentz.immstr.pasteur.fr/pdb/frozen_submission.php.

Additionally the frustratometer requires some external libraries to be installed in order to run.

You will need to install a library called BALL
http://www.ball-project.org/Downloads

And for BaLL to work you need to install:
http://www.ball-project.org/Downloads/index_html/Contrib/contrib-1.4.1.tar.gz
(if you use linux this is the link, see the link relative to BALL, for other OS or other linux versions)

You will also need to set an environmental variable in your .bashrc file like:
export BALL_DATA_PATH=/opt/BALL/BALL-1.4.0/data/
export BALL=${BALL}:/opt/BALL/BALL-1.4.0/build/include/BALL:/opt/BALL/BALL-1.4.0/build/BALL

changing the paths to the ones in your computer.
