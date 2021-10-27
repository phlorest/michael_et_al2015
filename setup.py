from setuptools import setup


setup(
    name='cldfbench_michael_et_al2015',
    py_modules=['cldfbench_michael_et_al2015'],
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'cldfbench.dataset': [
            'michael_et_al2015=cldfbench_michael_et_al2015:Dataset',
        ]
    },
    install_requires=[
        'phlorest',
        'ete3',
    ],
    extras_require={
        'test': [
            'pytest-cldf',
        ],
    },
)
