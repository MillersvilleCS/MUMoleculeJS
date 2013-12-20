(function(window) {
    'use strict';
    
    var SDFParser = function( ) {

    };

    SDFParser.parse = function(str) {

        var atoms = [];
        var bonds = [];

        var lines = str.split("\n");
        if (lines.length < 4)
            return;

        var atomCount = parseInt(lines[3].substr(0, 3));
        if (isNaN(atomCount) || atomCount <= 0)
            throw new IOException("Error loading SDF file, no atoms");

        var bondCount = parseInt(lines[3].substr(3, 3));
        if (lines.length < 4 + atomCount + bondCount)
            throw new IOException("Error loading SDF file, file too short");

        var offset = 4; //atoms start on line 4

        for (var i = 0; i < atomCount; i++) { //parse atoms
            var line = lines[offset];
            ++offset;

            atoms.push(parseAtom(line));
        }

        for (i = 0; i < bondCount; i++) { //parse bonds
            var line = lines[offset];
            ++offset;
            
            bonds.push(parseBond(line));
        }

        return {
            'ok': true,
            'atoms': atoms,
            'bonds': bonds
        };
    };

    function parseBond(line) {
        
        var from = parseInt(line.substr(0, 3));
        var to = parseInt(line.substr(3, 3));
        var count = parseInt(line.substr(6, 3));

        return {from: from - 1, to: to - 1, count: count};
    }

    function parseAtom(line) {
        
        var x = parseFloat(line.substr(0, 10));
        var y = parseFloat(line.substr(10, 10));
        var z = parseFloat(line.substr(20, 10));
        var e = trim(line.substr(30, 2)).toLowerCase();

        return {x: x, y: y, z: z, element: e, hflag: true};
    }

    function trim(text) {
        
        return text.replace(/^\s\s*/, '').replace(/\s\s*$/, '');
    }
    
    window.SDFParser = SDFParser;
})(window);
