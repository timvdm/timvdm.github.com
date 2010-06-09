goog.provide("jchemhub.view.Renderer");
goog.require("goog.structs.Map");


/**
 * Abstract Class to render a model object to a graphics object
 * 
 * @constructor
 * @param controller {goog.events.EventTarget} controller for this view
 * @param graphics {goog.graphics.AbstractGraphics} graphics to draw on.
 * @param opt_config {object} config to override defaults
 * @defaultConfig {object} object holding default values
 */
jchemhub.view.Renderer = function(controller, graphics, opt_config, defaultConfig) {
	this.controller = controller;
	this.graphics = graphics;

	this.config = new goog.structs.Map(defaultConfig);
	if (opt_config) {
		this.config.addAll(opt_config); // merge optional config into
		// defaults
	}
}

jchemhub.view.Renderer.prototype.render = goog.abstractMethod;
