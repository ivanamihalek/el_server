/*************************************************************************
    This code is from Dynamic Web Coding at dyn-web.com
    Copyright 2011 by Sharon Paine 
    See Terms of Use at www.dyn-web.com/business/terms.php
    regarding conditions under which you may use this code.
    This notice must be retained in the code as is!
*************************************************************************/

// change content onmouseover

// constructor
function dw_ContentSwap( opts ) { // className, displayId, restoreDefault, content
    this.cls = opts['className'] || 'showInfo'; 
    this.id = opts['displayId'] || 'infoDiv'; 
    this.restoreDefault = (typeof opts['restoreDefault'] === 'boolean')? opts['restoreDefault']: true;
    if ( this.restoreDefault ) {
        this.defaultInfo = document.getElementById( this.id ).innerHTML;
    }
    this.content = opts['content'];
    if ( this.content ) {
        this.init();
    } else {
        throw new Error('It seems you have not provided content to be displayed for dw_ContentSwap.');
    }
}

// for setting up events
dw_ContentSwap.prototype.init = function () {
    var els = dw_Util.getElementsByClassName( this.cls );
    var c, id, info;
    for (var i=0; els[i]; i++) {
        c = dw_Util.getNextClass( els[i].className, this.cls );
        
        if ( this.content[c] ) {
            id = this.id; info = this.content[c];
            dw_Event.add( els[i], 'mouseover', function( id, info ){
                return function (e) {
                    dw_ContentSwap.doSwap( id, info );
                }
            }( id, info ) ); // see Crockford js good parts pg 39
            
            if ( this.restoreDefault ) { // mouseout?
                info = this.defaultInfo; 
                dw_Event.add( els[i], 'mouseout', function( id, info ){
                    return function (e) {
                        var tgt = dw_Event.getTarget(e);
                        if ( dw_Util.mouseleave(e, tgt) ) {
                            dw_ContentSwap.doSwap( id, info );
                        }
                    }
                }( id, info ) );
            }
        }
    }
}

// class methods
dw_ContentSwap.doSwap = function (id, info) {
    document.getElementById(id).innerHTML = info;
}

// calls constructor, like rotator setup
dw_ContentSwap.setup = function () {
    for (var i=0; arguments[i]; i++) {
        new dw_ContentSwap( arguments[i] );
    }
}

/////////////////////////////////////////////////////////////////////
// Helper functions 

var dw_Util; if (!dw_Util) dw_Util = {};

// what className attached to what element type in what container element (default: document)
dw_Util.getElementsByClassName = function (sClass, sTag, oCont) {
    var result = [], list, i;
    var re = new RegExp("\\b" + sClass + "\\b", "i");
    oCont = oCont? oCont: document;
    if ( document.getElementsByTagName ) {
        if ( !sTag || sTag == "*" ) { // for ie5
            list = oCont.all? oCont.all: oCont.getElementsByTagName("*");
        } else {
            list = oCont.getElementsByTagName(sTag);
        }
        for (i=0; list[i]; i++) 
            if ( re.test( list[i].className ) ) result.push( list[i] );
    }
    return result;
}

// if actuatorClass (showTip, showInfo) in class name, next class would point to content 
dw_Util.getNextClass = function (cls, act) {
    if (!cls) return ''; var c = '';
        var classes = cls.split(/\s+/);
        for (var i=0; classes[i]; i++) {
            if ( classes[i] == act && classes[i + 1] ) {
                c = classes[i + 1];
                break;
            }
        }
    return c; // return next class name or ''
}

dw_Util.mouseleave = function (e, oNode) {
    e = dw_Event.DOMit(e);
    var toEl = e.relatedTarget? e.relatedTarget: e.toElement? e.toElement: null;
    if ( oNode != toEl && !dw_Util.contained(toEl, oNode) ) {
        return true;
    }
    return false;
}

dw_Util.contained = function (oNode, oCont) {
    if (!oNode) return null; // in case alt-tab away while hovering (prevent error)
    while ( (oNode = oNode.parentNode) ) if ( oNode == oCont ) return true;
    return false;
}

