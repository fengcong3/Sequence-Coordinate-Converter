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

setup(
    name='SequenceCoordinateConverter',
    version='0.1',
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
