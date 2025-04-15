from setuptools import setup, Extension, find_packages

# C Extension
genome2_processor_cext = Extension(
    'bri',  # name of the extension
    sources=[
        'scc/bri/bri_get.c',
        'scc/bri/bri_index.c',
        'scc/bri/py_bri_wrapper.c'  # include the wrapper
    ],  # source files
    extra_compile_args=[
        '-O3',
        '-std=c99',
        '-fsigned-char',
        '-D_FILE_OFFSET_BITS=64',
        '-g',
        '-fPIC'
    ],
    extra_link_args=[
        '-L/public/home/fengcong/anaconda2/lib/',
        '-fPIC',
        '-L/public/agis/chengshifeng_group/fengcong/WGRS/software/htslib-1.18/install_08284/lib/'
    ],
    include_dirs=[
        '/public/agis/chengshifeng_group/fengcong/WGRS/software/htslib-1.18/install_08284/include/',
        '/public/home/fengcong/anaconda2/envs/py3/include/python3.7m/'
    ],
    library_dirs=[
        '/public/agis/chengshifeng_group/fengcong/WGRS/software/htslib-1.18/install_08284/lib/'
    ],
    libraries=['pthread', 'z', 'm', 'lzma', 'hts']
)

# 如果有额外的包含路径，确保添加进去
ext_modules = [
    Extension('scc.bri.bri',
              sources=['scc/bri/py_bri_wrapper.c', 'scc/bri/bri_index.c', 'scc/bri/bri_get.c'],
              include_dirs=['/public/agis/chengshifeng_group/fengcong/WGRS/software/htslib-1.18/install_08284/include/'],
              libraries=['hts'],
              library_dirs=['/public/agis/chengshifeng_group/fengcong/WGRS/software/htslib-1.18/install_08284/lib/'],
              extra_compile_args=['-O3', '-std=c99', '-fsigned-char', 
                                 '-D_FILE_OFFSET_BITS=64', '-g', '-fPIC']
             )
]

setup(
    name='SequenceCoordinateConverter',
    version='0.3',  # 更新版本号
    packages=find_packages(),
    install_requires=[
        'pysam',
        'pandas'
    ],
    ext_modules=[genome2_processor_cext],
    entry_points={
        'console_scripts': [
            'SCC=scc.main:main',
        ],
    },
)
