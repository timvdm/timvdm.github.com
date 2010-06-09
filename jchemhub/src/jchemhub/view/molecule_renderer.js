goog.provide('jchemhub.view.MoleculeRenderer');
goog.require('jchemhub.controller.BondController');
goog.require('jchemhub.view.BondRenderer');
goog.require('jchemhub.view.BondRendererFactory');
goog.require('jchemhub.view.AtomRenderer');
goog.require('jchemhub.controller.AtomController');

/**
 * Class to render a molecule object to a graphics object
 * 
 * @constructor
 * @param graphics
 *            {goog.graphics.AbstractGraphics} graphics to draw on.
 * @extends {jchemhub.view.Renderer}
 */
jchemhub.view.MoleculeRenderer = function(controller, graphics, opt_config) {
	jchemhub.view.Renderer.call(this, controller, graphics, opt_config,
			jchemhub.view.MoleculeRenderer.defaultConfig);
	this.scale_factor = 1;
	this.bondController = new jchemhub.controller.BondController(controller);
	this.bondRendererFactory = new jchemhub.view.BondRendererFactory(
			this.bondController, graphics);
	this.atomController = new jchemhub.controller.AtomController(controller);
	this.atomRenderer = new jchemhub.view.AtomRenderer(this.atomController,
			graphics);
}
goog.inherits(jchemhub.view.MoleculeRenderer, jchemhub.view.Renderer);

jchemhub.view.MoleculeRenderer.prototype.render = function(molecule, trans,
		group) {
	if (!trans) {
		// if not part of a reaction, we need to create a transform
		trans = this.getTransform(molecule);
	}
	this.transform = trans;
	goog.array.forEach(molecule.bonds, function(bond) {
		this.bondRendererFactory.get(bond).render(bond, trans, undefined);
	}, this);
	goog.array.forEach(molecule.atoms, function(atom) {
		this.atomRenderer.render(atom, trans, this.atomController);
	}, this);
}

/**
 * 
 * @param {jchemhub.model.Molecule}
 *            reaction
 * @return {jchemhub.graphics.AffineTransform}
 */
jchemhub.view.MoleculeRenderer.prototype.getTransform = function(molecule) {
	var coords = goog.array.map(molecule.atoms, function(a) {
		return a.coord;
	})
	var m = Number(this.config.get("margin"));
	var box = goog.math.Box.boundingBox.apply(null, coords);
	var fromRect = goog.math.Rect.createFromBox(box.expand(m, m, m, m));

	var toSize = fromRect.getSize().scaleToFit(this.graphics.getSize());
	var scale = this.scale_factor * toSize.width / fromRect.getSize().width;

	var transform = new jchemhub.graphics.AffineTransform(scale, 0, 0, -scale,
			-fromRect.left * scale, -fromRect.top * scale);

	return transform;
};

/**
 * A default configuration for renderer
 */
jchemhub.view.MoleculeRenderer.defaultConfig = {
	margin : 4
};