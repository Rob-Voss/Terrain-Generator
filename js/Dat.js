/**
 *
 * @param object
 * @name datGUI
 * @constructor
 */
function datGUI (object) {
    var wGUI = new dat.GUI({resizable: true, autoPlace: true, name: 'Controls'}),
        sFolder = wGUI.addFolder('Settings');

    for (let go in object.settings) {
        if (object.settings.hasOwnProperty(go)) {
            let ty = typeof object.settings[go];
            if (ty === 'number') {
                sFolder.add(object.settings, go).listen().name(go).onFinishChange(function(){
                    object.drawMap();
                }).step(0.1);
            } else if (ty !== 'object') {
                sFolder.add(object.settings, go).listen().name(go).onFinishChange(function(){
                    object.drawMap();
                });
            } else {
               let folder = sFolder.addFolder(go.charAt(0).toUpperCase() + go.slice(1));
                for (var gp in object.settings[go]) {
                    if (object.settings[go].hasOwnProperty(gp)) {
                        let typ = typeof object.settings[go][gp];
                        if (typ !== 'object') {
                            folder.add(object.settings[go], gp).listen().name(gp).onFinishChange(function(){
                                object.drawMap();
                            });
                        }
                    }
                }
            }
        }
    }
    var fFolder = wGUI.addFolder('Functions');
    for (let fo in object.functions) {
        if (object.functions.hasOwnProperty(fo)) {
            let ty = typeof object.functions[fo];
            if (ty === 'function') {
                fFolder.add(object.functions, fo).listen().name(fo);
            } else {
                let folder = fFolder.addFolder(fo.charAt(0).toUpperCase() + fo.slice(1));
                for (var gp in object.functions[fo]) {
                    if (object.functions[fo].hasOwnProperty(gp)) {
                        let typ = typeof object.functions[fo][gp];
                        if (typ === 'function') {
                            folder.add(object.functions[fo], gp).listen().name(gp);
                        }
                    }
                }
            }
        }
    }
    wGUI.remember(object);
}
