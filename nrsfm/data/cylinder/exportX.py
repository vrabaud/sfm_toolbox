#!BPY
import Blender
from Blender import Mesh

filePath = '/users/u1/vrabaud/SSFM/data/cylinder/coord.txt'
out = open( filePath, 'w')

for iFrame in range(1,201):
	Blender.Set( 'curframe', iFrame )
	curve = Blender.Object.Get( 'Curve' )
	out.write( '%i Inf Inf\n' % ( iFrame ) )
	
	mesh = Mesh.New()
	mesh.getFromObject(curve.name)
	
	for vert in mesh.verts :
		out.write( '%f %f %f\n' % ( vert.co.x, vert.co.y, vert.co.z ) )

out.close()
