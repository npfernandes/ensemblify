ensemblify.config
=================

.. py:module:: ensemblify.config

.. autoapi-nested-parse::

   View and update Ensemblify's global configuration.



Attributes
----------

.. autoapisummary::

   ensemblify.config.GLOBAL_CONFIG


Functions
---------

.. autoapisummary::

   ensemblify.config.show_config
   ensemblify.config.update_config


Module Contents
---------------

.. py:data:: GLOBAL_CONFIG

   Ensemblify's global configuration dictionary.

.. py:function:: show_config() -> dict

   Show the current configuration dictionary.

   :returns:     The Ensemblify global configuration dictionary.
   :rtype: dict


.. py:function:: update_config(new_config: dict)

   Update the configuration dictionary with new values.

   Any keys in the new configuration that exist in the
   current config are updated.

   :param new_config: Dictionary with configuration parameters to update.
   :type new_config: :py:class:`dict`


