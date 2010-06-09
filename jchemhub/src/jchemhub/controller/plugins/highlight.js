
goog.provide('jchemhub.controller.plugins.Highlight');
goog.require('jchemhub.controller.Plugin');
goog.require('goog.functions');
goog.require('goog.debug.Logger');


/**
 * simple Plugin for highlighting bonds and atoms
 *
 * @constructor
 * @extends {jchemhub.controller.Plugin}
 */
jchemhub.controller.plugins.Highlight = function() {
  jchemhub.controller.Plugin.call(this);
};
goog.inherits(jchemhub.controller.plugins.Highlight, jchemhub.controller.Plugin);


/** @inheritDoc */
jchemhub.controller.plugins.Highlight.prototype.getTrogClassId =
    goog.functions.constant('jchemhub.controller.plugins.Highlight');

/**
 * Logging object.
 * 
 * @type {goog.debug.Logger}
 * @protected
 */
jchemhub.controller.plugins.Highlight.prototype.logger = goog.debug.Logger
		.getLogger('jchemhub.controller.plugins.Highlight');


jchemhub.controller.plugins.Highlight.prototype.handleAtomMouseOver = function(e) {
	if (!e.currentTarget.atomHighlightGroup) {
		e.currentTarget.atomHighlightGroup = this.highlightAtom(e.atom);
	} else {
		e.currentTarget.atomHighlightGroup = this.highlightAtom(e.atom, e.currentTarget.atomHighlightGroup);
	}

};

jchemhub.controller.plugins.Highlight.prototype.handleAtomMouseOut = function(e) {

	if (e.currentTarget.atomHighlightGroup) {
		e.currentTarget.atomHighlightGroup.clear();
		e.currentTarget.atomHighlightGroup = undefined;
	}

};

jchemhub.controller.plugins.Highlight.prototype.handleBondMouseOver = function(e) {
	this.logger.info("handleBondMouseOver");
	if (!e.currentTarget.bondHighlightGroup) {
		e.currentTarget.bondHighlightGroup = this.highlightBond(e.bond);
	} else {
		e.currentTarget.bondHighlightGroup = this.highlightBond(e.bond, e.currentTarget.bondHighlightGroup);
	}

};
jchemhub.controller.plugins.Highlight.prototype.handleBondMouseOut = function(e) {

	if (e.currentTarget.bondHighlightGroup) {
		e.currentTarget.bondHighlightGroup.clear();
		e.currentTarget.bondHighlightGroup = undefined;
	}

};

jchemhub.controller.plugins.Highlight.prototype.highlightBond = function(bond, opt_group){
	return this.editorObject.reactionRenderer.moleculeRenderer.bondRendererFactory.get(bond).highlightOn(bond, opt_group);
};

jchemhub.controller.plugins.Highlight.prototype.highlightAtom = function(atom, opt_group){
	return this.editorObject.reactionRenderer.moleculeRenderer.atomRenderer.highlightOn(atom, opt_group);
};
