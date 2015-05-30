# -*- mode: python -*-
a = Analysis(['Gui.py'],
             pathex=['C:\\Users\\Nick\\Documents\\GitHub\\Voxtel-Dither'],
             hiddenimports=[],
             hookspath=None,
	     excludes=['PySide'],
             runtime_hooks=None)
for d in a.datas:
    if 'pyconfig' in d[0]: 
        a.datas.remove(d)
        break
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='PhasePlateDesigner.exe',
          debug=False,
          strip=None,
          upx=True,
          console=False )
