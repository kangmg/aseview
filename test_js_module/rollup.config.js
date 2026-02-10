import terser from '@rollup/plugin-terser';

const banner = `/**
 * ASEView.js v0.1.0
 * JavaScript library for molecular visualization
 * https://github.com/kangmg/aseview
 * (c) 2024 MIT License
 */`;

export default [
    // UMD build (for browsers via CDN)
    {
        input: 'src/index.js',
        output: {
            file: 'dist/aseview.umd.js',
            format: 'umd',
            name: 'ASEView',
            banner,
            globals: {
                three: 'THREE'
            }
        },
        external: ['three']
    },
    // Minified UMD build
    {
        input: 'src/index.js',
        output: {
            file: 'dist/aseview.umd.min.js',
            format: 'umd',
            name: 'ASEView',
            banner,
            globals: {
                three: 'THREE'
            }
        },
        external: ['three'],
        plugins: [terser()]
    },
    // ES module build
    {
        input: 'src/index.js',
        output: {
            file: 'dist/aseview.esm.js',
            format: 'es',
            banner
        },
        external: ['three']
    }
];
