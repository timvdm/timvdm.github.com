goog.provide('jchemhub.controller.MoleculeController');
goog.require('goog.events.EventTarget');
goog.require('goog.debug.Logger');

/** 
 * @constructor 
 * @extends {goog.events.EventTarget} 
 */ 
jchemhub.controller.MoleculeController = function(parentController) { 
  goog.events.EventTarget.call(this);
  this.setParentEventTarget(parentController);

}; 
goog.inherits(jchemhub.controller.MoleculeController, goog.events.EventTarget); 

jchemhub.controller.MoleculeController.prototype.handleMouseOver = function(Molecule, e){
	this.dispatchEvent(jchemhub.controller.MoleculeController.EventType.MOUSEOVER);
};

jchemhub.controller.MoleculeController.prototype.handleMouseOut = function(Molecule, e){
	this.dispatchEvent(jchemhub.controller.MoleculeController.EventType.MOUSEOUT);
};
/** @enum {string} */ 
jchemhub.controller.MoleculeController.EventType = { 
  MOUSEOVER: 'molecule_mouseover',
  MOUSEOUT: 'molecule_mouseout'
}; 

/**
 * Logging object.
 * 
 * @type {goog.debug.Logger}
 * @protected
 */
jchemhub.controller.MoleculeController.prototype.logger = goog.debug.Logger
		.getLogger('jchemhub.controller.MoleculeController');
